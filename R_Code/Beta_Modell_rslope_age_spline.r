library(nimble)
library(coda)
library(splines2)
## if directly run use this lines first to get the data setup
data <- readRDS("..\\data\\AIdataset_normalized.Rds")
data$id <- as.numeric(data$id)
summary(data$ai)
data$ai <- data$ai + 0.001

#we want to regress: (sofa ~ age + day +(1|id), data) in a bayesian way
summary(data$age)


#has 6 values 
#spline for  
subjects <- length(unique(data$id))
model.inits <- list(
  list(beta0 = 0,   beta_age=rep(0, 6), beta_day = 0,    phi = 5,  sigma_b0 = 1.0, sigma_b1 = 1, b0 = rep(0.0,times = subjects), b1 = rep(0.0,times = subjects)),
  list(beta0 = 0.5, beta_age=rep(1, 6), beta_day = -0.1, phi = 10, sigma_b0 = 1.5, sigma_b1 = 1, b0 = rep(0.2,times = subjects), b1 = rep(0.0,times = subjects)),
  list(beta0 = 1,   beta_age=rep(1, 6), beta_day = 0.1,  phi = 7,  sigma_b0 = 0.5, sigma_b1 = 1, b0 = rep(0.3,times = subjects), b1 = rep(0.0,times = subjects))
)
parameters = c("beta0", "beta_age", "beta_day","sigma_b0","sigma_b1","phi","Deviance","ppo")
model.constants <- list( N = length(data$sofa), x1 = data$age,
                         id = data$id,          x2 = data$day,  
                         Nsubj = length(unique(data$id)))

#instead of the normal x include b-splines
bsMat_age <- bSpline(data$age, knots = c(-20,0,20),degree = 3)

model.data <- list(y = data$ai, X_age = bsMat_age)

model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dbeta(a[i], b[i])
    a[i] <- mu[i] * phi
    b[i] <-(1 - mu[i]) * phi
    #choice of a l logit link for this parameter
    #beta1 age beta_day day
    logit(mu[i]) <- beta0 + inprod(beta_age[], X_age[i,]) +
      beta_day *x2[i]+ b0[id[i]]+ x2[i] *b1[id[i]]
    ppo[i] <- dbeta(y[i],a[i], b[i])
    D[i] <- -2*log(ppo[i])
  }
  Deviance <- sum(D[1:N])
  #priors
  phi ~ dunif(0,10000)
  sigma_b0 ~ dgamma(0.001, 0.001)
  sigma_b1 ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,sd = 1000)
  beta_day ~ dnorm(0,sd = 1000)
  for(j in 1:6){
    beta_age[j] ~ dnorm(0, tau = 5)
  }
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,sd = sigma_b0)
    b1[i] ~ dnorm(0,sd = sigma_b1)
  }
})


t_0 <- Sys.time()
ai_spline_age <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 11000,
  nburnin = 5000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC
t_1 <- Sys.time()
time_linear <- t_1 - t_0
time_linear
saveRDS(ai_spline_age, file = "..\\data\\mcmc_res\\ai_spline_age_rslope.Rds")
