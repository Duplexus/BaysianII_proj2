library(nimble)
library(coda)
## if directly run use this lines first to get the data setup
data <- readRDS("..\\data\\AIdataset_normalized.Rds")
data$id <- as.numeric(data$id)
summary(data$ai)
data$ai <- data$ai + 0.001
#we want to regress: (sofa ~ age + day +(1|id), data) in a bayesian way

# model.inits <- list(list(phi=2,sigma_b0 =2, beta0=1, beta1 = 1,beta2 = 1,b0 = c(rep(1,times = length(unique(data$id))))),
#                     list(phi=20,sigma_b0 = 1, beta0=10, beta1 = 10,beta2 = 10, b0 =  rnorm(length(unique(data$id)),0,30)),
#                     list(phi=15, sigma_b0 =2, beta0=20, beta1 = -10,beta2 = -15, b0 =  runif(length(unique(data$id)),0,10))
# )

model.inits <- list(
  list(beta0 = 0, beta1 = 0, beta2 = 0, phi = 30, sigma_b0 = 1),
  list(beta0 = 0.5, beta1 = -0.1, beta2 = -0.1, phi = 10, sigma_b0 = 1.5),
  list(beta0 = 1, beta1 = 0.1, beta2 = 0.1, phi = 20, sigma_b0 = 0.5)
)
parameters = c("beta0", "beta1", "beta2","sigma_b0","phi","y_pred")

model.constants <- list( N = length(data$sofa), x1 = data$age,
                         id = data$id, x2 = data$day,  Nsubj = length(unique(data$id)))
model.data <- list(y = data$ai)

model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dbeta(a[i], b[i])
    y_pred[i] ~ dbeta(a[i], b[i])
    a[i] <- mu[i] * phi
    b[i] <-(1 - mu[i]) * phi
    #choice of a l logit link for this parameter
    logit(mu[i]) <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
  }
  #priors
  phi ~ dunif(0,10000)
  sigma_b0 ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,sd = 1000)
  beta1 ~ dnorm(0,sd = 1000)
  beta2 ~ dnorm(0,sd = 1000)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,sd = sigma_b0)
  }
})

beta_rand <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 5000,
  nburnin = 2000, thin = 1, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC

data$fail

