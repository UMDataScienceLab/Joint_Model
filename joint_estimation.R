
#' Simulation

require(Matrix) ## create diagonal matrix
require(nloptr) ## optimizer
require(numDeriv) ## gradient
require(MASS) ## basic package
library(survival) ## survival model
library(mnormt) ## multivaraite normal
library(nlme) ## mixed effect model
library(glmmML) ## Gauss-Hermite quadrature
library(AlgDesign)  ## factorial design
library(statmod)  ## Gauss-Legendre quadrature
require(ggplot2) # plots
library(pracma) ## Trace function
library(rootSolve) ## root solver
library(pROC) ## roc curve
require(boot) ## glm
library(e1071) ## svm

n_iter <- 1 ## number of iteration (experiment)

ROC_Matrix_true <- matrix(0,nrow = n_iter, ncol = 6) ## true prob
ROC_Matrix_hat <- matrix(0,nrow = n_iter, ncol = 6) ## estimated prob from MGP

for(iter in 1:n_iter){
  
  source("new_rejectSam.R") ## generate time-to-event data
  
  ######################################################################################
  ## Simulation Study 1
  ## 1 latent function X
  ## input is in one-dimension space (i.e. time)
  ######################################################################################
  
  ########################################### System Input ##############################################
  
  sigma <- sigma2
  n <- length(table(bdata$battid))+1 # Number of output  
  len <- unique(as.vector(table(bdata$battid))) # length of each output 
  m <- part_obs+1 # observations from last unit
  
  ########################################### Generate Data #############################################
  
  train_id <- unique(bdata$battid)
  trains <- lapply(train_id, function(i) bdata[bdata$battid==i,]$obstime)
  tests <- 0:part_obs
  trainy <- lapply(train_id, function(i) bdata[bdata$battid==i,]$r)
  testy <- test_data$r[(tests+1)]
  sum_pi <- sum(lengths(trains))+length(testy)
  
  ########################################## Plot ###############################################
  
  for(i in 1:(n-1)) {
    plot(trains[[i]],trainy[[i]],xlim=c(0,60),ylim=c(0,20),type="l",lwd=0.1,xlab=NA,ylab=NA,cex.axis=1.2)
    par(new=T)
  }
  par(new=T)
  plot(tests,testy,xlim=c(0,60),ylim=c(0,20),type="b",pch=2,col="blue",lwd=2,xlab=NA,ylab=NA,cex.axis=1.2)

  hada=bdiag(matrix(1,length(trains[[1]]),length(trains[[1]])))
  for(i in 2:(n-1))
    hada=bdiag(hada, matrix(1,length(trains[[i]]),length(trains[[i]])) )
  hada=bdiag(hada, length(tests), length(tests))
  index=as.matrix(which(hada!=0,arr.ind = T));
  
  ##################### Covariance Functions #####################################################
  ## 1) Latent
  cuu=function(a,b,L) {d=outer(a,b,`-`);I=outer(a,b,`==`);res=L[1]^2; exp(-0.5*d^2/res) +I*L[2*n+3]^2}
  
  ## 2) Function with Latent
  cu=function(a,b,L) {d=outer(a,b,`-`); res=L[1]^2+L[3]^2; L[2]*sqrt(L[1]^2/res)*exp(-0.5*d^2/res)}
  
  ## 3) Function with Function
  cf=function(a,b,L) {
    d=outer(a,b,`-`);I=outer(a,b,`==`)
    res=L[1]^2+2*L[3]^2
    L[2]^2 * sqrt(L[1]^2/res)
  }
  
  cf_cox=function(a,b,L) {
    d=outer(a,b,`-`);I=outer(a,b,`==`)
    res=L[1]^2+2*L[3]^2
    L[2]^2 * sqrt(L[1]^2/res) * exp(-0.5*d^2/res)
  }
  
  ##################### Covariance (Faster through Sparse matrix) ####################################
  cfu=function(L){
    
    temp <- do.call(  rbind, lapply(1:(n-1), function(i)cu(trains[[i]],z,L[c(1,i+1,n+1+i)]))  )
    return(rbind(temp, cu(tests,z,L[c(1,n+1,n+1+n)]) ))
    
  }
  
  cff=function(L)
  {
    # only diagonal elements of cff are non-zero since we only consider trace
    
    ff <- do.call(sum, lapply(1:(n-1), function(i) length(trains[[i]])*cf(trains[[i]], trains[[i]], L[c(1,i+1,n+1+i)])) )
    ff <- ff + length(tests) * cf(tests, tests, L[c(1,n+1,n+1+n)]) 
    
    return(ff)
  }
  
  integrand <- function(t, L, kk, r_cuu, m_vec, S_vec){
    as.numeric(t>floor(min(V)))*exp(L[2*n+4]+abs(L[2*n+5])*(t-min(V)) )*exp(L[2*n+6]*max(bdata[bdata$battid==train_id[kk],]$x)+
          L[2*n+7]*(cu(t,z,L[c(1,kk+1,n+1+kk)])%*%solve(r_cuu,m_vec,tol=1e-32)+0.5*diag( cf_cox(t,t,L[c(1,kk+1,n+1+kk)]) -
                      cu(t,z,L[c(1,kk+1,n+1+kk)])%*%solve(r_cuu, diag(1,length(z))-S_vec%*%solve(r_cuu, tol=1e-32), tol=1e-32)%*%t(cu(t,z,L[c(1,kk+1,n+1+kk)])) )  ))
  }
  
  integrand_test <- function(t, L, kk, r_cuu, m_vec, S_vec){
    
    as.numeric(t> floor(min(V)))*exp(L[2*n+4]+abs(L[2*n+5])*(t-min(V)) )*exp(L[2*n+6]*max(test_data$x)+
          L[2*n+7]*(cu(t,z,L[c(1,kk+1,n+1+kk)])%*%solve(r_cuu,m_vec,tol=1e-32)+0.5*diag( cf_cox(t,t,L[c(1,kk+1,n+1+kk)]) -
                      cu(t,z,L[c(1,kk+1,n+1+kk)])%*%solve(r_cuu, diag(1,length(z))-S_vec%*%solve(r_cuu, tol=1e-32), tol=1e-32)%*%t(cu(t,z,L[c(1,kk+1,n+1+kk)])))  ))
    }

  
  
  cox_1 <- function(L, r_cuu, m_vec, S_vec){
    
    ## L: (2n+4):(2n+7), b, phi, gamma, beta
   
    cox <- do.call(sum, lapply(1:(n-1), function(i) 
      max(bdata[bdata$battid==train_id[i],]$status)*(L[2*n+4]+abs(L[2*n+5])*(V[train_id[i]]-min(V) )+L[2*n+6]*max(bdata[bdata$battid==train_id[i],]$x)+
                                                       L[2*n+7]*cu(V[train_id[i]],z,L[c(1,i+1,n+1+i)])%*%solve(r_cuu,m_vec,tol=1e-32)) 
      ))
    cox <- cox + 0 
    
    cox2 <- do.call(sum, lapply(1:(n-1), function(kk) integrate(f=integrand, lower = 0, upper = c(V[train_id[kk]]), L, kk, r_cuu, m_vec, 
                                                                S_vec, subdivisions=1000, rel.tol =.Machine$double.eps^0.25*10, stop.on.error = F)$value ))
    cox2 <- cox2 + integrate(f=integrand_test, lower = 0, upper = tail(tests, 1), L, n, r_cuu, m_vec, S_vec, subdivisions=1000, 
                             rel.tol = .Machine$double.eps^0.25*10, stop.on.error = F)$value
    return(cox-cox2)
  }
  
  
  ##################### Liklihood ##################################################################
  
  y=c(unlist(trainy), testy) ## input signals
  z=seq(0,40,length.out=10) ## pseudo-input
  
  logL=function(H,fn)
  {
    
    ## log-marginal likelihood
    
    r_cuu <- cuu(z,z,H); r_cfu <- cfu(H); r_cff <- cff(H)
    fuf <- r_cfu%*%solve(r_cuu,t(r_cfu), tol=2e-32)
    D <- r_cff - sum(diag(fuf)) 
    B <- fuf + diag(H[2*n+2]^2, dim(fuf)[1])
    
    m_vec <- H[2*n+2]^(-2)*r_cuu%*%solve(r_cuu+H[2*n+2]^(-2)*t(r_cfu)%*%r_cfu, t(r_cfu), tol=1e-32)%*%y 
    S_vec <-r_cuu%*%solve(r_cuu+H[2*n+2]^(-2)*t(r_cfu)%*%r_cfu, r_cuu, tol=1e-32)
    
    deter <- det(B)
    if(deter>0) {
      a <- 0.5*(log(deter)+t(y)%*%solve(B,y,tol=2e-32)+log(2*pi)*length(y)+H[2*n+2]^(-2)*D) - cox_1(H, r_cuu, m_vec, S_vec)
    } else {
      ch <- chol(B, pivot=TRUE, tol=2e-32)
      logdeter <- 2*(sum(log(diag(ch))))
      a <- 0.5*(logdeter+t(y)%*%solve(B,y,tol=2e-32)+log(2*pi)*length(y)+H[2*n+2]^(-2)*D) - cox_1(H, r_cuu, m_vec, S_vec)
    }
    return(as.numeric(a))
  }
  
  logL_grad=function(H,fn) {return(nl.grad(H,fn))} ## gradient
  
  ## Method of asymptotes : check all options in list using / nloptr.print.options() 
  x0 <- c(1,rep(1,2*n),0.1,0.5, 0.1, 0.1, 0.1, 0.5)
  opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 500, print_level = 3)
  one <- nloptr(x0=x0,eval_f = logL,eval_grad_f = logL_grad,opts = opts,fn = logL, lb=c(rep(-Inf, 2*n+3), -Inf, -Inf, -Inf, -Inf ),
                ub=c(rep(Inf, 2*n+3), Inf, Inf, Inf, Inf ))
  H0 <- one$solution ## optimal solution 
  
  ####################### Results #################################################################
  
  #1) Input Cov
  r_cuu=cuu(z,z,H0);r_cfu=cfu(H0);r_cff=cff(H0)
  fuf=r_cfu%*%solve(r_cuu,t(r_cfu), tol=2e-32)
  D=sparseMatrix(i=index[,1],j=index[,2],x=(r_cff-fuf)[index])
  B=fuf + diag(H0[2*n+2]^2, dim(fuf)[1])
  
  #2) Mean 
  r1_main=solve(B,y); r2_main=solve(r_cuu,t(r_cfu), tol=2e-32)
  
  #3) calculate fitted mean for all individuals
  temp_matrix <- NULL
  for(i in 1:length(train_id)) temp_matrix <- rbind(temp_matrix, as.matrix(cu(trains[[i]],z,H0[c(1,i+1,i+n+1)])%*%r2_main%*%r1_main) )
  bdata$rfit = temp_matrix; rm(temp_matrix)
  
  risk_time_set <-unique(  sort(c(c(1: max(V)), V))  )
  H0t_est_cumsum <- cumsum(exp(H0[2*n+4]+H0[2*n+5]*(risk_time_set-min(V) )) )
  H0t_sm <- smooth.spline(risk_time_set, H0t_est_cumsum, nknots=4, spar = 0.4)
  gam0 = H0[2*n+6]
  beta11 = H0[2*n+7]
  
  #4) check fitted MGP
  for(i in 1:length(train_id)){
    ypred = cu(trains[[i]],z,H0[c(1,i+1,i+n+1)])%*%r2_main%*%r1_main
    plot(trains[[i]],ypred,xlim=c(0,30),ylim=c(0,15),"l",lty=1,lwd=1,cex.axis=1.2,ann=FALSE)
    par(new=T)
    plot(trains[[i]],trainy[[i]],xlim=c(0,30),ylim=c(0,15),"l",lty=2,lwd=1,cex.axis=1.2,ann=FALSE)
  }
  
  #5) event probability estimation
  Bp = Bmat[testid, ]
  xp=test_data$x[1]
  iit = 1
  tstar = part_obs
  PAR = list(H0t_sm=H0t_sm, beta10=beta11, gam0=gam0)
  
  Pr1 = 1 - St_cond(tstar+8,tstar,x=xp,Bp)   #probability for the battery failure within one year from tstar
  Pr2 = 1 - St_cond(tstar+10,tstar,x=xp,Bp)
  Pr3 = 1 - St_cond(tstar+12,tstar,x=xp,Bp)
  Pr4 = 1 - St_cond(tstar+15,tstar,x=xp,Bp)
  Pr5 = 1 - St_cond(tstar+20,tstar,x=xp,Bp)
  
  ROC_Matrix_true[iter,] <-c(Pr1, Pr2, Pr3, Pr4, Pr5, max(test_data$obstime)-part_obs)
  
  Pr1_EST = 1 - St_cond_EST(tstar+8,tstar,x=xp,PAR)
  Pr2_EST = 1 - St_cond_EST(tstar+10,tstar,x=xp,PAR)
  Pr3_EST = 1 - St_cond_EST(tstar+12,tstar,x=xp,PAR)
  Pr4_EST = 1 - St_cond_EST(tstar+15,tstar,x=xp,PAR)
  Pr5_EST = 1 - St_cond_EST(tstar+20,tstar,x=xp,PAR)
  
  ROC_Matrix_hat[iter,] <- c(Pr1_EST, Pr2_EST, Pr3_EST, Pr4_EST, Pr5_EST, max(test_data$obstime)-part_obs)
  
  #6) mean and 95CI
  r3_main=r2_main%*%solve(B,t(r2_main))
  ypred3=cu(test_data$obstime,z,H0[c(1,n+1,2*n+1)])%*%r2_main%*%r1_main
  yvar3=abs(diag( cf(test_data$obstime,test_data$obstime,H0[c(1,n+1,2*n+1)])-cu(test_data$obstime,z,H0)%*%r3_main%*%t(cu(test_data$obstime,z,H0[c(1,n+1,2*n+1)])) )  )
  
  plot(test_data$obstime,ypred3,xlim=c(0,30),ylim=c(0,15),"l",lty=1,lwd=1,cex.axis=1.2,ann=FALSE)
  par(new=T)
  plot(test_data$obstime,test_data$r,xlim=c(0,30),ylim=c(0,15),"l",lty=2,lwd=1,cex.axis=1.2,ann=FALSE,col="green")
  par(new=T)
  plot(test_data$obstime,ypred3+3*sqrt(yvar3),lwd=1,xlim=c(0,30),ylim=c(0,15),"l",lty=2,cex.axis=1.2,ann=FALSE)
  par(new=T)
  plot(test_data$obstime,ypred3-3*sqrt(yvar3),lwd=1,xlim=c(0,30),ylim=c(0,15),"l",lty=2,cex.axis=1.2,ann=FALSE)

}


