#ex 03 first try:
library("coda")
library("survival")
library(nimble)
library(lattice)
library(dplyr)
source("helpfunctions.R")
dataman <- readRDS("../AIdataset.RDS")
pred_sofa <- readRDS("../y_pred4_sofa.Rds")
pred_ai <- readRDS("../pred_ai.Rds")
data_e03 <- readRDS("Prep_data_for_ex3.rds")
data_e03$day <- as.numeric(dataman$day)
data_e03$age <- as.numeric(dataman$age)
data_e03$sofa_pred <- as.numeric(pred_sofa)
data_e03$ai_pred <- as.numeric(pred_ai)
data_e03$id <- as.numeric(data_e03$id)
data_e03$fail <- as.numeric(data_e03$fail)

data_e03$fail2 <- ifelse(data_e03$fail == 2, 0, data_e03$fail)
cens <- ifelse(data_e03$fail == 0, 30, 0)
is.censored <- ifelse(data_e03$fail2 == 0,1,0)
data_e03$t <- ifelse(is.censored == 1,NA,data_e03$day)

model.inits <- list(list(k=2, beta0=0.1, beta1 = 5,beta2 = 10 ),
                    list(k=1, beta0=20, beta1 = 20,beta2 = 20 ),
                    list(k=0.5, beta0=2, beta1 = -1,beta2 = -1 )
)

model.constants <- list( N = length(data_e03$fail), Nsubj=length(unique(data_e03$id)),id=data_e03$id)

data_e03$fail3 <- ifelse(data_e03$fail == 1, data_e03$icustaymv, 1000)
model.data <- list(y = data_e03$t, cens=cens, x1 = data_e03$sofa_pred,
                   x2 = data_e03$ai_pred)



# MODEL SPECIFICATION 
#now the model specification is in line with the estimation in survival
model.function <- nimbleCode({
  for (i in 1:N){
    censored[i] ~ dinterval(y[i],cens[i])
    y[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]] #b0 random intercept wrt patient
    # Distribution of future observed counts for plate i
    #predict[i]  ~ dweib(k, invlambda[i]) 
  }
  #priors
  scale <- 1/k
  k ~ dunif(0,1000)
  beta0 ~ dnorm(0,0.00001)
  beta1 ~ dnorm(0,0.00001)
  beta2 ~ dnorm(0,0.00001)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,1000)
  for(z in 1:Nsubj){
    b0[z] ~ dnorm(0,tau_b0)
  }
})


#Monitored Variables
parameters <-c("beta0", "beta1", "beta2","scale")

runjags.options(method = "rjparallel")
#Set Up Model

weibull_e03 <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 100000,
  nburnin = 20000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, 
  WAIC = T, summary=T) 
summarise_default(weibull_e03$samples)




#library(survival)
#data_e03$fail3 <- ifelse(data_e03$fail == 1, data_e03$icustaymv, 1000)
#(survi.aft <- survreg(formula = Surv(fail3) ~ sofa_pred+ ai_pred, data = data_e03, dist = "weibull"))
##survival time chagnes by this factor
#exp(coef(survi.aft))

#(shapeParameter <- 1 / survi.aft$scale)
##hazard rate changes by this factor
#exp(-1 * shapeParameter * coef(survi.aft))



# till here some Weibull interpretation 