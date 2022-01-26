library(nimble)
library(coda)
library(splines2)
## if directly run use this lines first to get the data setup
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#try(setwd(".\\data"))
source("..\\R_Code\\helpfunctions.r")
data <- readRDS("..\\data\\AIdataset_normalized.Rds")
data$id <- as.numeric(data$id)
summary(data$ai)
summary(data$day)
data$ai <- data$ai + 0.001

#we want to regress: (sofa ~ age + day +(1|id), data) in a bayesian way
summary(data$age)


#has 6 values 
#spline for  
groups_id <- length(unique(data$id))
model.inits <- list(
  list(lambda = 2, beta0 = 0.0, beta_age=0,beta_day=rep(0, 20),phi = 30, sigma_b0 = 1.0,b0 = rep(0,groups_id)),
  list(lambda = 2,beta0 = 0.5, beta_age=1,beta_day=rep(1, 20),phi = 10, sigma_b0 = 1.5,b0 = rep(0,groups_id)),
  list(lambda = 2,beta0 = 1.0, beta_age=1,beta_day=rep(1, 20),phi = 20, sigma_b0 = 0.5,b0 = rep(0,groups_id))
)
parameters = c("beta0", "beta_age", "beta_day","sigma_b0","phi","lambda")
library(splines2)

#instead of the normal x include b-splines
bsMat_day <- bSpline(data$day,degree = 3,df = 20)
penalization_mat <- makeQ(2,ncol(bsMat_day), epsilon = 0.1)
model.constants <- list( N = length(data$sofa),
                         id = data$id, 
                         K = ncol(bsMat_day),
                         Q = penalization_mat,
                         mean_penmat = rep(0,ncol(bsMat_day)),
                         Nsubj = length(unique(data$id)))

model.data <- list(y = data$ai, X_age = data$age, X_day = bsMat_day)

eigen(penalization_mat)$values

model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dbeta(a[i], b[i])
    a[i] <- mu[i] * phi
    b[i] <-(1 - mu[i]) * phi
    #choice of a l logit link for this parameter
    #beta1 age beta2 day
    logit(mu[i]) <-beta0+beta_age*X_age[i]+inprod(beta_day[], X_day[i,])+ b0[id[i]]
  }
  #priors
  phi ~ dunif(0,10000)
  sigma_b0 ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,sd = 1000)
  beta_age ~ dnorm(0,sd = 1000)
  lambda ~ dgamma(0.01, 0.01)
  A[1:K,1:K] <- lambda * Q[1:K,1:K]
  beta_day[1:K] ~ dmnorm(mean_penmat[1:K],A[1:K,1:K])
  for ( j in 1:Nsubj){
      b0[j] ~ dnorm(0,sd = sigma_b0)
  }
})
    
t_0 <- Sys.time()
ai_pbspline_day <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 10000,
  nburnin = 5000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, # get an coda object instead of plain values
  WAIC = T ) # get the WAIC
t_1 <- Sys.time()
time_linear <- t_1 - t_0
time_linear
saveRDS(ai_pbspline_day, file = "..\\data\\mcmc_res\\day_pen_spline.Rds")

