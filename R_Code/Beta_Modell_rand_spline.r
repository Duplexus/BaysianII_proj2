library(nimble)
library(coda)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

obs_per <- as.data.frame(data %>%  group_by(id) %>%  summarise(lange = length(id)))["lange"]
mean(unlist(as.vector(obs_per)))
#quantile(unlist(as.vector(obs_per)),c((0:10)*0.1))

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

model.inits <- list(
  list(beta0 = 0, beta_age = 0, beta_day = 0, phi = 30, sigma_b0 = 1),
  list(beta0 = 0.5, beta_age = -0.1, beta_day = -0.1, phi = 10, sigma_b0 = 1.5),
  list(beta0 = 1, beta_age = 0.1, beta_day = 0.1, phi = 20, sigma_b0 = 0.5)
)
parameters = c("beta0", "beta_age", "beta_day","sigma_b0","phi","Deviance","ppo","b0","y_pred")

model.constants <- list( N = length(data1$sofa), x1 = data1$age,
                         id = data1$id, x2 = data1$day,  Nsubj = length(unique(data1$id)))
model.data <- list(y = data1$ai)

model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dbeta(a[i], b[i])
    #just take y_pred[i] if the results are good, (they are good looks like)
    y_pred[i] ~ dbeta(a[i], b[i])
    a[i] <- mu[i] * phi
    b[i] <-(1 - mu[i]) * phi
    #choice of a l logit link for this parameter
    logit(mu[i]) <- beta0 + beta_age *x1[i] + beta_day *x2[i]+ b0[id[i]]
    ppo[i] <- dbeta(y[i],a[i], b[i])
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
  }
})

beta_rand <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 25000,
  nburnin = 15000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC

saveRDS(beta_rand, file = "..\\data\\mcmc_res\\beta_rand.rds")
plot(10:100,dgamma(10:100,0.001, 0.001),type = "l")


sigma_b0 <- runif(1000,0,10000)
tau_b0 <- 1/sigma_b0
hist(1/(tau_b0))

