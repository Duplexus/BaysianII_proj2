library(nimble)
library(coda)
library(splines2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## if directly run use this lines first to get the data setup
data2 <- readRDS("..\\data\\AIdataset_normalized.Rds")
data2$id <- as.numeric(data2$id)
summary(data2$ai)
data2$ai <- data2$ai + 0.001

#we want to regress: (sofa ~ age + day +(1|id), data) in a bayesian way
summary(data2$age)


#has 6 values 
#spline for  
groups_id <- length(unique(data2$id))
model.inits <- list(
  list(beta0 = 0.0, beta_age=0,beta_day=rep(0, 6),phi = 30, sigma_b0 = 1.0,b0 = rep(0,groups_id)),
  list(beta0 = 0.5, beta_age=1,beta_day=rep(1, 6),phi = 10, sigma_b0 = 1.5,b0 = rep(0,groups_id)),
  list(beta0 = 1.0, beta_age=1,beta_day=rep(1, 6),phi = 20, sigma_b0 = 0.5,b0 = rep(0,groups_id))
)
parameters = c("beta0", "beta_age", "beta_day","sigma_b0","phi","Deviance","ppo","b0")
model.constants <- list( N = length(data2$sofa),
                         id = data2$id, 
                         Nsubj = length(unique(data2$id)))
library(splines2)

#instead of the normal x include b-splines

dd_quant <- quantile(data2$day,c(0.25,0.5,0.75))
bsMat_day <- bSpline(data2$day, knots = c(dd_quant[1],dd_quant[2],dd_quant[3]),degree = 3)
model.data <- list(y = data2$ai, X_age = data2$age, X_day = bsMat_day)

model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dbeta(a[i], b[i])
    a[i] <- mu[i] * phi
    b[i] <-(1 - mu[i]) * phi
    #choice of a l logit link for this parameter
    #beta1 age beta2 day
    logit(mu[i]) <-beta0+beta_age*X_age[i]+inprod(beta_day[], X_day[i,])+ b0[id[i]]
    ppo[i] <- dbeta(y[i],a[i], b[i])
    D[i] <- -2*log(ppo[i])
  }
  Deviance <- sum(D[1:N])
  #priors
  phi ~ dunif(0,10000)
  sigma_b0 ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,sd = 1000)
  beta_age ~ dnorm(0,sd = 1000)
  for(j in 1:6){  
    beta_day[j] ~ dnorm(0, tau = 5)
  }
  for ( j in 1:Nsubj){
      b0[j] ~ dnorm(0,sd = sigma_b0)
  }
})
    
ai_spline_day <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 20000,
  nburnin = 10000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC
saveRDS(ai_spline_day, file = "..\\data\\mcmc_res\\ai_spline_day.Rds")

rm(data2)
