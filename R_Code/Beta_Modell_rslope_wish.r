library(nimble)
library(coda)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## if directly run use this lines first to get the data setup
data1 <- readRDS("..\\data\\AIdataset_normalized.Rds")
data1$id <- as.numeric(data1$id)
summary(data1$ai)
data1$ai <- data1$ai + 0.001
subjects <- length(unique(data1$id))

model.inits <- list(
  list(beta0 = 0,   beta_age=0, beta_day = 0,    phi = 5, b0 = rep(0.0,times = subjects), b1 = rep(0.0,times = subjects)),
  list(beta0 = 0.5, beta_age=1, beta_day = -0.1, phi = 10, b0 = rep(0.2,times = subjects), b1 = rep(0.0,times = subjects)),
  list(beta0 = 1,   beta_age=1, beta_day = 0.1,  phi = 7,  b0 = rep(0.3,times = subjects), b1 = rep(0.0,times = subjects))
)
parameters = c("beta0", "beta_age", "beta_day","phi","Deviance","ppo",
               "b1","b0","sigma2_11","sigma2_12","sigma2_22")

R <- diag(c(1,2))
model.constants <- list( N = length(data1$sofa), x1 = data1$age, R = R,
                         id = data1$id, x2 = data1$day,  Nsubj = length(unique(data1$id)))
model.data <- list(y = data1$ai)

model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dbeta(a[i], b[i])
    #just take y_pred[i] if the results are good, (they are good looks like)
    a[i] <- mu[i] * phi
    b[i] <-(1 - mu[i]) * phi
    #choice of a l logit link for this parameter
    logit(mu[i]) <- beta0 + beta_age *x1[i] + beta_day *x2[i]+ b0[id[i]]+ beta_day *b1[id[i]]
    ppo[i] <- dbeta(y[i],a[i], b[i])
    D[i] <- -2*log(ppo[i])
  }
  #priors
  Deviance <- sum(D[1:N])
  phi ~ dunif(0,10000)
  beta0 ~ dnorm(0,sd = 1000)
  beta_age ~ dnorm(0,sd = 1000)
  beta_day ~ dnorm(0,sd = 1000)
  sigma2_b[1:2,1:2] ~ dinvwish(S = R[1:2,1:2], df = 2)
  sigma2_11 <- sigma2_b[1,1]
  sigma2_12<- sigma2_b[1,2]
  sigma2_22 <- sigma2_b[2,2]
  #it becomes more vague as smaller the diag and as so small the df
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,sd = sigma2_b[1,1])
    b1[i] ~ dnorm(0,sd = sigma2_b[2,2])
  }
})

beta_rand <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 10000,
  nburnin = 5000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC

saveRDS(beta_rand, file = "..\\data\\mcmc_res\\beta_rand_rslope_wish2.rds")

beta_rand$WAIC$WAIC

