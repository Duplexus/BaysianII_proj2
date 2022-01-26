#Excercise 4

library(nimble)
library(coda)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## if directly run use this lines first to get the data setup
data1 <- readRDS("..\\data\\AIdataset_normalized.Rds")
data1$id <- as.numeric(data1$id)
summary(data1$ai)
data1$ai <- data1$ai + 0.001
#we want to regress: (sofa ~ age + day +(1|id), data) in a bayesian way

# model.inits <- list(list(phi=2,sigma_b0 =2, beta0=1, beta_age = 1,beta_day = 1,b0 = c(rep(1,times = length(unique(data1$id))))),
#                     list(phi=20,sigma_b0 = 1, beta0=10, beta_age = 10,beta_day = 10, b0 =  rnorm(length(unique(data1$id)),0,30)),
#                     list(phi=15, sigma_b0 =2, beta0=20, beta_age = -10,beta_day = -15, b0 =  runif(length(unique(data1$id)),0,10))
# )

model.inits1 <- list(
  list(beta0 = 0, beta_age = 0, beta_day = 0, phi = 30, sigma_b0 = 1),
  list(beta0 = 0.5, beta_age = -0.1, beta_day = -0.1, phi = 10, sigma_b0 = 1.5),
  list(beta0 = 1, beta_age = 0.1, beta_day = 0.1, phi = 20, sigma_b0 = 0.5)
)
parameters1 = c("beta0", "beta_age", "beta_day","sigma_b0","phi","Deviance","ppo","b0","y1_pred")

model.constants1 <- list( N = length(data1$sofa),
                         id = data1$id,  Nsubj = length(unique(data1$id)))
model.data1 <- list(y1 = data1$ai, x2 = data1$day,x1 = data1$age)

model.function <- nimbleCode({
  #Model for AI
  for (i in 1:N){
    y1[i] ~ dbeta(a[i], b[i])
    #just take ya_pred[i] if the results are good, (they are good looks like)
    y1_pred[i] ~ dbeta(a[i], b[i])
    a[i] <- mu[i] * phi
    b[i] <-(1 - mu[i]) * phi
    #choice of a l logit link for this parameter
    logit(mu[i]) <- beta0 + beta_age *x1[i] + beta_day *x2[i]+ b0[id[i]]
    ppo[i] <- dbeta(y1[i],a[i], b[i])
    D[i] <- -2*log(ppo[i])
  }
  #priors
  Deviance <- sum(D[1:N])
  phi ~ dunif(0,10000)
  sigma_b0 ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,sd = 1000)
  beta_age ~ dnorm(0,sd = 1000)
  beta_day ~ dnorm(0,sd = 1000)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,sd = sigma_b0)
    
    
  # Model for SOFA
    
    
    
  # Model for combination
    
    
    
  }
})

beta_rand2 <- nimbleMCMC(
  code = model.function, constants = model.constants1,
  data = model.data1, inits = model.inits1,
  monitors = parameters1, nchains = 3, niter = 25000,
  nburnin = 15000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC

#saveRDS(beta_rand, file = "..\\data\\mcmc_res\\beta_rand.rds")

