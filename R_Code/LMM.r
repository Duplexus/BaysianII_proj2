library(nimble)
library(coda)
## if directly run use this lines first to get the data setup
data <- readRDS("..\\data\\AIdataset_normalized.Rds")
data$id <- as.numeric(data$id)
#we want to regress: (sofa ~ age + day +(1|id), data) in a bayesian way

model.inits <- list(list(tau=2,sigma2_b0 =2, beta0=1, beta1 = 1,beta2 = 1,b0 = c(rep(1,times = length(unique(data$id))))),
                    list(tau=20,sigma2_b0 = 1, beta0=10, beta1 = 10,beta2 = 10, b0 =  rnorm(length(unique(data$id)),0,30)),
                    list(tau=15, sigma2_b0 =2, beta0=20, beta1 = -10,beta2 = -15, b0 =  runif(length(unique(data$id)),0,10))
)
parameters = c("sigma2", "beta0", "beta1", "beta2","sigma2_b0","tau","y_pred")

model.constants <- list( N = length(data$sofa), x1 = data$age,
                    id = data$id, x2 = data$day,  Nsubj = length(unique(data$id)))
model.data <- list(y = data$sofa)

model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dnorm(mu[i], tau)
    y_pred[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
  }
  #priors
  tau ~ dgamma(0.001,0.001)
  sigma2 <- (1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
})

sofa_lmm <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 5000,
  nburnin = 2000, thin = 1, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC

data$fail

