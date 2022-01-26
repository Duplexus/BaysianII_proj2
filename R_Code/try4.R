#ex 03 first try:
library("coda")
library("survival")
library(nimble)
library(lattice)
library(dplyr)
library("statmod")
dataman <- readRDS("AIdataset.RDS")

n <- length(unique(dataman$id))
M <- table(dataman$id)

#surv cens times
dataman$fail2 <- ifelse(dataman$fail == 2, 0, dataman$fail)
cens <- ifelse(dataman$fail == 0, 30, 0)
is.censored <- ifelse(dataman$fail2 == 0,1,0)
dataman$t <- ifelse(is.censored == 1,NA,dataman$day)
Time <- matrix(n)
treat <- matrix(n)

#set $t time as censoring time ie 30 if fail=0 or max(day) if fail=1
for(i in 1:n){
  tmp_t=ifelse(dataman[dataman$id==i,]$fail2 == 0,max(dataman[dataman$id==i,"day"]),30)
  dataman[dataman$id==i,]$t <- tmp_t
  Time[i] <- unique(tmp_t)
  treat[i] <- unique(dataman[dataman$id==i,]$icustaymv)

  }

death <- dataman$fail2

#longitudutunal info in matrix
time <- matrix(NA,n,max(M))
di.data <- matrix(NA,n,max(M))
count <- 1
for(i in 1:n){
  di.data[i, 1:M[i]] <- dataman$di[count:(M[i]+count-1)]
  time[i,1:M[i]] <- dataman$day[count:(M[i]+count-1)]
  count <- count+M[i]
}

#split longitudnal into fixed XL and random ZL
XL <- array(1,dim=c(n,max(M),3)) #fix
XL[, , 2] <- time
XL[, , 3] <- treat
ZL <- array(1,dim=c(n,max(M),2)) #random effect
ZL[,,2] <- time
XS <- model.matrix(~treat)
#quadrature
glq <- gauss.quad(15,kind="legendre")
xk <- glq$nodes  #nodes
wk <- glq$weights #weights
K <- length(xk) #k points



model.function <- nimbleCode({
  for(i in 1:n){
    #longit
    for(j in 1:M[i]){
      di.data[i,j] ~ dnorm(mu[i,j],tau)
      mu[i,j] <- inprod(betaL[],XL[i,j,1:3])+inprod(b[i,1:2],ZL[i,j,1:2])
    }
    #surv and cens
    for(j in 1:K){
      haz[i,j] <- alpha*pow(Time[i]/2*(xk[j]+1),alpha-1) *
        exp(inprod(betaS[],XS[i,1:2])+gamma*(b[i,1]+b[i,2]*(Time[i]/2*(xk[j]+1))))
    }
    #logsurv w gauss-legendre quad
    logSurv[i] <- -Time[i]/2 * inprod(wk,haz[i,1:15])
    
    #surv loglikehliodo
    phi[i] <- 100000 - death[i]*log(haz[i,K]) - logSurv[i]
    zeros[i]~dpois(phi[i])
    
    #random effects
    b[i,1:Nb] ~ dmnorm(mub[],Omega[,])
  }
  #priors
  for(l in 1:NbetasS){
    betaS[l]~dnorm(0,0.0001)
  }
  gamma~dnorm(0,0.0001)
  alpha~dunif(0,100)
  for(l in 1:NbetasL){
    betaL[l]~dnorm(0,0.0001)
  }
  tau <- pow(sigma,-2)
  sigma~dunif(0,100)
  Omega[1:Nb,1:Nb]~dwish(V[,],Nb)
  Sigma[1:Nb,1:Nb]<- inverse(Omega[,])
  
  lambda <- exp(betaS[1])
})



##elemnt list
d.jags <- list(n = n, M = M, Time = Time, XS = XS, 
               di.data = di.data, XL = XL, ZL = ZL, 
               death = death, mub = rep(0, 2), V = diag(1, 2), 
               Nb = 2, zeros = rep(0, n), 
               NbetasL = dim(XL)[3], NbetasS = ncol(XS), 
               K = length(xk), xk = xk, wk = wk)

i.jags <- function(){
  list(betaS=rnorm(ncol(XS)), gamma=rnorm(1),alpha=runif(1),
       betaL=rnorm(dim(XL)[3]),sigma=runif(1),Omega=diag(runif(2)))
}

p.jags <- c("betaS","gamma","alpha","lambda","betaL","sigma","Sigma","b")



library("rjags")

#Monitored Variables

#runjags.options(method = "rjparallel")
#Set Up Model

consts = list(Nb=2,NbetasS=ncol(XS),NbetasL=dim(XL)[3])
try4 <- nimbleMCMC(
  code = model.function, constants=consts,
  data = d.jags, inits = i.jags,
  monitors = p.jags, nchains = 3, niter = 20000,
  nburnin = 3000, thin = 3, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, 
  WAIC = T, summary=T) 

try4 <- jags.model(data=d.jags,file="try4mod.txt",inits=i.jags,n.chains=3)
try4.res <- coda.samples(try4,variable.names=p.jags,n.iter=10000,
                         thin=3)

#summary(weibull_e03$samples)
result <- as.mcmc(do.call(rbind,try4.res))
Sigma2.11 <- result[,1]; Sigma2.12 <- result[,2]; Sigma2.22 <- result[,4]
alpha <- result[,5]; b1 <- result[,6:(n+5)]; b2 <- result[, (n+6):(2*n+5)]
betaL1 <- result[,(2*n+6)]; betaL2 <- result[,(2*n+7)]; betaL3 <- result[,(2*n+8)]
betaS2 <- result[,(2*n+10)]; gamma <- result[,(2*n+11)]
lambda <- result[,(2*n+12)]; sigma <- result[,(2*n+13)]


#library(survival)
#data_e03$fail3 <- ifelse(data_e03$fail == 1, data_e03$icustaymv, 1000)
#(survi.aft <- survreg(formula = Surv(fail3) ~ sofa_pred+ ai_pred, data = data_e03, dist = "weibull"))
##survival time chagnes by this factor
#exp(coef(survi.aft))

#(shapeParameter <- 1 / survi.aft$scale)
##hazard rate changes by this factor
#exp(-1 * shapeParameter * coef(survi.aft))



# till here some Weibull interpretation 