
##----------------------------------------------------------------------------------
##      The integrand G(z) in the survival prediction
##----------------------------------------------------------------------------------

## for GP
testFunction_pred <- function(t0){
  xstar <- t0
  r_cuu=cuu(z,z,H0);r_cfu=cfu(H0);r_cff=cff(H0)
  fuf=r_cfu%*%solve(r_cuu,t(r_cfu), tol=2e-32)
  D=sparseMatrix(i=index[,1],j=index[,2],x=(r_cff-fuf)[index])
  B=fuf + diag(H0[2*n+2]^2, dim(fuf)[1])
  r1_main=solve(B,y); r2_main=solve(r_cuu,t(r_cfu), tol=2e-32)
  ypred = cu(xstar,z,H0[c(1,n+1,2*n+1)])%*%r2_main%*%r1_main
  return(ypred)
}


ht_sm = function(t,H0t_sm,beta1,gam0,x)  # survival function
{ 
  as.numeric(predict(H0t_sm,t,deriv=1)$y*exp(gam0*x+ beta1*(testFunction_pred(t)))) 
}

G = function (H0t_sm, tstar, t, beta1,gam0,x)
{
  tp = ht_sm(seq(tstar,t,length.out=1000),H0t_sm,beta1,gam0,x) #G in formula (11)
  return(exp(-sum(tp*(t-tstar)/999)))
}


##----------------------------------------------------------------------------------
##      The function below keeps both the return value and warnings generated from a function
##----------------------------------------------------------------------------------  

keepWarnings <- function(expr) 
{
   localWarnings <- list()
   value <- withCallingHandlers(expr,
   warning = function(w) {
   localWarnings[[length(localWarnings)+1]] <<- w
   invokeRestart("muffleWarning")
   })
   list(value=value, warnings=localWarnings)
}



