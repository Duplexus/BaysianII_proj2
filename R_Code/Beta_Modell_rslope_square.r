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
subjects <- length(unique(data1$id))

model.inits <- list(
  list(beta0 = 0,   beta_age=0, beta_day = 0,    phi = 5, sigma_b0 = 1.0, sigma_b1 = 1, sigma_b2 = 1,b0 = rep(0.0,times = subjects), b1 = rep(0.0,times = subjects), b2 = rep(0.0,times = subjects)),
  list(beta0 = 0.5, beta_age=1, beta_day = -0.1, phi = 10, sigma_b0 = 1.5, sigma_b1 = 1,sigma_b2 = 1,b0 = rep(0.2,times = subjects), b1 = rep(0.0,times = subjects), b2 = rep(0.0,times = subjects)),
  list(beta0 = 1,   beta_age=1, beta_day = 0.1,  phi = 7, sigma_b0 = 0.5, sigma_b1 = 1, sigma_b2 = 1,b0 = rep(0.3,times = subjects), b1 = rep(0.0,times = subjects), b2 = rep(0.0,times = subjects))
)
parameters = c("beta0", "beta_age", "beta_day","sigma_b0","sigma_b1","sigma_b2","phi","Deviance","ppo","b0")

model.constants <- list( N = length(data1$sofa), x1 = data1$age,
                         id = data1$id, x2 = data1$day,  Nsubj = length(unique(data1$id)))
model.data <- list(y = data1$ai)

model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dbeta(a[i], b[i])
    a[i] <- mu[i] * phi
    b[i] <-(1 - mu[i]) * phi
    logit(mu[i]) <- beta0 + beta_age *x1[i] + beta_day *x2[i]+ b0[id[i]]+ 
      x2[i] *b1[id[i]]+ pow(x2[i],2) *b2[id[i]]
    ppo[i] <- dbeta(y[i],a[i], b[i])
    D[i] <- -2*log(ppo[i])
  }
  #priors
  Deviance <- sum(D[1:N])
  phi ~ dunif(0,10000)
  sigma_b0 ~ dgamma(0.001, 0.001)
  sigma_b1 ~ dgamma(0.001, 0.001)
  sigma_b2 ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,sd = 1000)
  beta_age ~ dnorm(0,sd = 1000)
  beta_day ~ dnorm(0,sd = 1000)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,sd = sigma_b0)
    b1[i] ~ dnorm(0,sd = sigma_b1)
    b2[i] ~ dnorm(0,sd = sigma_b2)
  }
})

beta_rand <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 14000,
  nburnin = 7000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC

saveRDS(beta_rand, file = "..\\data\\mcmc_res\\beta_rand_rslope_square.rds")

