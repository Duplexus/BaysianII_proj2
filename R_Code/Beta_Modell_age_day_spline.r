library(nimble)
library(coda)
library(splines2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## if directly run use this lines first to get the data setup
data3 <- readRDS("..\\data\\AIdataset_normalized.Rds")
data3$id <- as.numeric(data3$id)
summary(data3$ai)
data3$ai <- data3$ai + 0.001

#we want to regress: (sofa ~ age + day +(1|id), data) in a bayesian way
summary(data3$age)


#has 6 values 
#spline for  
groups_id <- length(unique(data3$id))
model.inits <- list(
  list(beta0 = 0.0, beta_age=rep(0, 6),beta_day=rep(0, 6),phi = 30, sigma_b0 = 1.0,b0 = rep(0,groups_id)),
  list(beta0 = 0.5, beta_age=rep(1, 6),beta_day=rep(1, 6),phi = 10, sigma_b0 = 1.5,b0 = rep(0,groups_id)),
  list(beta0 = 1.0, beta_age=rep(1, 6),beta_day=rep(1, 6),phi = 20, sigma_b0 = 0.5,b0 = rep(0,groups_id))
)
parameters = c("beta0", "beta_age", "beta_day","sigma_b0","phi","ppo","Deviance")
model.constants <- list( N = length(data3$sofa),
                         id = data3$id, 
                         Nsubj = length(unique(data3$id)))
library(splines2)

#instead of the normal x include b-splines
#da_quant <- quantile(data3$age,c(0.25,0.5,0.75))
bsMat_age <- bSpline(data3$age, knots = c(-20,0,20),degree = 3)
dd_quant <- quantile(data3$day,c(0.25,0.5,0.75))
bsMat_day <- bSpline(data3$day, knots = c(dd_quant[1],dd_quant[2],dd_quant[3]),degree = 3)
model.data <- list(y = data3$ai, X_age = bsMat_age, X_day = bsMat_day)

model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dbeta(a[i], b[i])
    a[i] <- mu[i] * phi
    b[i] <-(1 - mu[i]) * phi
    #choice of a l logit link for this parameter
    #beta1 age beta2 day
    logit(mu[i]) <- beta0 + inprod(beta_age[], X_age[i,]) +
      inprod(beta_day[], X_day[i,])+ b0[id[i]]
    ppo[i] <- dbeta(y[i],a[i], b[i])
    D[i] <- -2*log(ppo[i])
  }
  Deviance <- sum(D[1:N])
  #priors
  phi ~ dunif(0,10000)
  sigma_b0 ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,sd = 1000)
  for(j in 1:6){
    beta_age[j] ~ dnorm(0, tau = 5)
    beta_day[j] ~ dnorm(0, tau = 5)
  }
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,sd = sigma_b0)
  }
})


t_0 <- Sys.time()
ai_spline_age_day <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 20000,
  nburnin = 10000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC
t_1 <- Sys.time()
time_linear <- t_1 - t_0
time_linear

saveRDS(ai_spline_age_day, file = "..\\data\\mcmc_res\\ai_spline_age_day.rds")
