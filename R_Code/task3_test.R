#ex 03 first try:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("coda")
library("survival")
library(nimble)
library(lattice)
library(dplyr)
source("helpfunctions.R")
dataman <- readRDS("../data/AIdataset_normalized.Rds")
pred_sofa <- readRDS("../data/y_pred4_sofa.Rds")
pred_ai <- readRDS("../data/pred_ai.Rds")
data_e03 <- dataman

data_e03$pred_ai <- pred_ai 
data_e03$pred_sofa <- pred_sofa 

data_e03 <- data_e03 %>%  group_by(id) %>%  mutate(maxdays = max(day))
data_e03 <- data_e03[data_e03$maxdays == data_e03$day,]

data_e03$fail2 <- ifelse(data_e03$fail == 2, 0, data_e03$fail)
data_e03$true_surv <- ifelse(data_e03$fail2  == 1,data_e03$icustaymv,30)

data_e03 <- data_e03 %>% arrange(true_surv,fail2)


model.inits <- list(list(k=2, beta0=0.1, beta_sofa = 5,beta_ai = 10 ),
                    list(k=1, beta0=20, beta_sofa = 20,beta_ai = 20 ),
                    list(k=0.5, beta0=2, beta_sofa = -1,beta_ai = -1 )
)

censored_times <- ifelse(data_e03$true_surv == 30,NA, data_e03$true_surv )
censored <- ifelse(data_e03$true_surv == 30, 1, 0)


model.constants <- list( N = length(data_e03$fail),id=as.numeric(data_e03$id))

model.data <- list(y = censored_times, censored=censored, x1 = as.numeric(data_e03$pred_sofa),
                   x2 = as.numeric(data_e03$pred_ai))
# model.data <- list(y = censored_times, censored=censored, x1 = as.numeric(data_e03$sofa),
#                    x2 = as.numeric(data_e03$ai))



# MODEL SPECIFICATION 
#now the model specification is in line with the estimation in survival
model.function <- nimbleCode({
  # for (i in 1:N1){
  #   y[i] ~ dweib(k, invlambda[i])
  #   invlambda[i] <- pow(t[i], k)
  #   t[i] <- exp(-h[i])
  #   h[i] <- beta0 + beta_sofa *x1[i] + beta_ai *x2[i]
  # }
  for (i in (0+1):N){
    censored[i] ~ dinterval(y[i],30)
    y[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta_sofa *x1[i] + beta_ai *x2[i]
  }
  #priors
  scale <- 1/k
  k ~ dunif(0,1000)
  beta0 ~ dnorm(0,0.00001)
  beta_sofa ~ dnorm(0,0.00001)
  beta_ai ~ dnorm(0,0.00001)
})



#Monitored Variables
parameters <-c("beta0", "beta_sofa", "beta_ai","scale")

weibull_e03 <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 10000,
  nburnin = 5000, thin = 10, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, 
  WAIC = T, summary=T) 
summarise_default(weibull_e03$samples)




#library(survival)
a <- (survi.aft <- survreg(formula = Surv(icustaymv) ~ pred_sofa+ pred_ai, data = data_e03, dist = "weibull"))
b <- (survi.aft <- survreg(formula = Surv(icustaymv) ~ sofa+ ai, data = data_e03, dist = "weibull"))

summary(b)
summary(a)
#survival time chagnes by this factor
# exp(coef(survi.aft))
# 
# (shapeParameter <- 1 / survi.aft$scale)
# #hazard rate changes by this factor
# exp(-1 * shapeParameter * coef(survi.aft))

# till here some Weibull interpretation