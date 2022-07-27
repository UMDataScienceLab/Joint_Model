

source("R_functions_new_new.R") # rejection sampling function; keepWarnings()

## Data Parameters
Ns=3 ## Number of Units
Niter=1 ## Number of Iterations

## Baseline Hazard Weibul
lambda = 0.001                   
alpha = 1.05

## Cox Parameters
gam = 0 
beta1 = 0.5 

## Degradation Signal Parameters (quadratic form)
coe <- seq(3,30, length.out = Ns)/1000
MUB = matrix(0, nrow = Ns, ncol = 3)
for(i in 1:nrow(MUB)){
  MUB[i, ] = c(2.5, 0.01, coe[i])
}
SIGMAB = matrix(c(0.2,-4E-4,-8E-5,-4E-4,3E-6,3E-7,-8E-5,3E-7,1E-7),3, 3, byrow=TRUE)
sigma2 = 0.5

## survival functions
h0t = function(t) lambda*alpha*(t^(alpha-1))   # Baseline hazard
H0t = function(t) lambda*(t^alpha)
S0t = function(t) exp(-lambda*(t^alpha))   # baseline survival

Gt = function(t,B) B[1]+B[2]*t+B[3]*t^2 # For generating resistance signal
ht_true = function(t,x,B) {h0t(t)*exp(gam*x + beta1*(Gt(t,B)))}
St_true_tp = function(t,x,B) exp(-integrate(f=ht_true, lower=1E-3, upper=t, x,B, subdivisions=1000)$value)
St_true = Vectorize(St_true_tp, "t")
ft_true = function(t,x,B) St_true(t,x,B)*ht_true(t,x,B) #survival event density function

## true survival function
St_cond = function(t,tstar,x,Bp) St_true(t,x,Bp)/St_true(tstar,x,Bp)

## define the function for estimated conditional survival function (Gauss Quad)
St_cond_EST_tp = function(t, tstar, x, PAR) 
{
  Itp = G(PAR$H0t_sm, tstar, t, PAR$beta10, PAR$gam0, x)
  return(Itp)
}
St_cond_EST = Vectorize(St_cond_EST_tp, "t")   # vectorised function

## rejection sampling
RejectionSampling = function(n,knownDensity,A,B,minRgDensity,maxRgDensity) ## Generate Death Times
{
  f = match.fun(knownDensity)  #now f is the density to sample from
  fvalues = f(seq(from=minRgDensity,to=maxRgDensity,length.out=100),A,B)
  if (sum(fvalues*(maxRgDensity-minRgDensity)/99)<0.95)
  {
    return(NA)
    stop
  }
  maxDensityValue = max(fvalues)*1.1 #the maximum value of knownDensity(x) 
  RN = NULL
  for(i in 1:n)
  {
    OK = 0
    while(OK<1)
    {
      T = runif(1,min = minRgDensity, max = maxRgDensity )
      U = runif(1,min = 0, max = 1)
      if(U*maxDensityValue <= f(T,A,B))
      {
        OK = 1
        RN = c(RN,T)
      }
    }
  }
  return(RN)
}

## Database Generation

generation=function()
{
  iit = 0
  while (1)
  {
    iit = iit+1
    cat(" \n \n")
    cat("STARTING ITERATION",iit,"...\n")
    
    Bmat <<- matrix(0, nrow = Ns, ncol = 3)
    for(i in 1:Ns){
      Bmat[i, ] <<- rmnorm(1, MUB[i, ], SIGMAB) 
    }
    
    X = rbinom(Ns, 1, 0.5)  # random indicators of x (1-> manufacturer A; 0-> B)
    T = c()  # generated death times
    for (i in 1:Ns)  T[i] = RejectionSampling(1,ft_true,X[i],Bmat[i,],minRgDensity=0.001,maxRgDensity=120)
    V <<- T
    
    NA_sampling_index <- which(is.na(V)) #resample if it is NA
    while(length(NA_sampling_index)!=0){
      for (i in NA_sampling_index) 
        V[i] <<- RejectionSampling(1,ft_true,X[i],Bmat[i,],minRgDensity=0.001,maxRgDensity=120)
      NA_sampling_index <- which(is.na(V))
    }
    
    #### Generate Censoring Data
    Cid = ceiling(runif(ceiling(Ns*0.05),0,Ns))  # id's of the 5% censored batteries
    for (i in 1:ceiling(Ns*0.05))
    {
      V[Cid[i]] <<- runif(1,min=0,max=T[Cid[i]])  # censoring time
    }
    
    status = rep(1,Ns)
    status[Cid] = 0
    
    tr = sum(ceiling(V))   # total number of rows in the data
    ER = rnorm(tr, mean = 0, sd = sqrt(sigma2))  # random measurement errors
    rNA = rep(NA,tr)
    # data structure: x-baseline covariate, r-resistance measurement, status-died(1)/censored(0),obstime-time instance r is made, eventtime-died/censored time
    # (start,stop,died,rini,rfit,rinifit) is only for initial parameter value estimation of the survival part, not for joint modeling
    bdata = data.frame(battid = rNA, x = rNA, r = rNA, status = rNA, obstime = rNA, eventtime = rNA, start = rNA, stop = rNA, died = rNA, rini = rNA, rfit = rNA, rinifit = rNA, nerr=rNA, nerrini=rNA)
    
    for (i in 1:Ns)
    {
      ir = ceiling(V)[i]  # number of rows for the ith battery
      #print(i)
      ts = sum(ceiling(V)[0:(i-1)])  # total number of rows in bdata for the 1~(i-1)th batteries
      bdata[(ts+1):(ts+ir),1] = rep(i,ir)    # battery id
      bdata[(ts+1):(ts+ir),2] = rep(X[i],ir)    # indicator of x (1-> manufacturer A; 0-> B)
      bdata[(ts+1):(ts+ir),13] = Gt(0:floor(V[i]),Bmat[i,]) # Signal without noise
      bdata[(ts+1):(ts+ir),3] = bdata[(ts+1):(ts+ir),13]+ER[(ts+1):(ts+ir)]   # observed resistance
      bdata[(ts+1):(ts+ir),4] = rep(status[i],ir)    # died(1)/censored(0)
      bdata[(ts+1):(ts+ir),5] = 0:floor(V[i])    # observation time points (assumption: every month)
      bdata[(ts+1):(ts+ir),6] = rep(V[i],ir)   # death/censored time point
      # elements below are for bdata_surv part (for only for initial parameter value estimation of the survival part)
      start = 0:(ir-1)
      stop = 1:ir
      stop[length(stop)] = V[i]
      died = c(rep(0,ir-1),status[i])
      bdata[(ts+1):(ts+ir),7] = start
      bdata[(ts+1):(ts+ir),8] = stop
      bdata[(ts+1):(ts+ir),9] = died
      bdata[(ts+1):(ts+ir),10] = rep(bdata[(ts+1),3],ir) # initial resistance
      bdata[(ts+1):(ts+ir),14] = rep(bdata[(ts+1),13],ir) # initial resistance without error
      
    }
    
    cat("     data generation completed...\n")
    if (iit==Niter) {cat(" \n","N = ",Niter,"is reached\n");break}
  }
  return(list(bdata=bdata,Bmat=Bmat,V=V))
}

check_data_quality <- function(){
  
  flag = 0
  counter = 0
  while(flag==0){
    
    counter <- counter +1
    data=generation();
    bdata=data$bdata;V=data$V
    if(length(unique(bdata$battid))>=2) ##geneally more than 2 units
      flag=1
  }
  
  print("good data!")
  print(counter)
  print("times")
  return(bdata)
  
}

bdata <- check_data_quality() ##get a good quality 
testid <- sample(unique(bdata$battid), size=1)
test_data <- bdata[bdata$battid==testid, ]
h <- ggplot(data=test_data, aes(x=obstime, y=r)) 
h + geom_line(aes(group = battid, colour = battid)) 
percent <- 0.6 #can vary from 0 to 1
Rr_true <- matrix(test_data$r, ncol=1)
part_obs <- floor(percent*length(Rr_true)) 
Rr <- as.matrix(Rr_true[1:(part_obs+1)], part_obs + 1, 1)
bdata <- bdata[bdata$battid!=testid, ]
h <- ggplot(data=bdata, aes(x=obstime, y=r)) 
h + geom_line(aes(group = battid, colour = battid)) 
RH <- matrix(bdata$r, ncol = 1) #historical data
bdata_surv=Surv(bdata$start,bdata$stop,bdata$died)
