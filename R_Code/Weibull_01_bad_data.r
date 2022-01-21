#ex 03 first try:
library("coda")
library("survival")
data_e03 <- readRDS("..\\data\\Prep_data_for_ex3.rds")
data_e03$fail2 <- ifelse(data$fail == 2, 0, data$fail)

model.inits <- list(list(k=3, beta0=1, beta1 = 1,beta2 = 1 ),
                    list(k=1, beta0=10, beta1 = 10,beta2 = 10 ),
                    list(k=0.5, beta0=20, beta1 = -10,beta2 = -15 )
)

model.constants <- list( N = length(data_e03$fail), x1 = data_e03$sofa_pred,
                         x2 = data_e03$ai_pred)

data_e03$fail3 <- ifelse(data_e03$fail == 1, data_e03$icustaymv, 1000)
model.data <- list(y = data_e03$fail3 )



# MODEL SPECIFICATION 
#now the model specification is in line with the estimation in survival
model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    # Distribution of future observed counts for plate i
    #predict[i]  ~ dweib(k, invlambda[i]) 
  }
  #priors
  scale <- 1/k
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
})


#Monitored Variables
parameters <-c("beta0", "beta1", "beta2","scale")

runjags.options(method = "rjparallel")
#Set Up Model

weibull_e03 <- nimbleMCMC(
  code = model.function, constants = model.constants,
  data = model.data, inits = model.inits,
  monitors = parameters, nchains = 3, niter = 5000,
  nburnin = 2000, thin = 1, setSeed = c(1,2,3),
  samplesAsCodaMCMC = T, samples = T, 
  WAIC = T ) 
summary(weibull_e03$samples)




library(survival)
data_e03$fail3 <- ifelse(data_e03$fail == 1, data_e03$icustaymv, 1000)
(survi.aft <- survreg(formula = Surv(fail3) ~ sofa_pred+ ai_pred, data = data_e03, dist = "weibull"))
#survival time chagnes by this factor
exp(coef(survi.aft))

(shapeParameter <- 1 / survi.aft$scale)
#hazard rate changes by this factor
exp(-1 * shapeParameter * coef(survi.aft))



# till here some Weibull interpretation 



