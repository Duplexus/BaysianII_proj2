setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("runjags")
library("coda")
library("rjags")
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
#attributes(Grub$grubsize) <- NULL


#DEFINE INTITIAL VALUES
model.inits <- list(list(tau=2, beta0=1, beta1 = 1,beta2 = 1 ),
                    list(tau=20, beta0=10, beta1 = 10,beta2 = 10 ),
                    list(tau=15, beta0=20, beta1 = -10,beta2 = -15 )
)

#Monitored Variables
parameters <-c("beta0", "beta1", "beta2", "sigma2","tau")



#what happens if big grub size is forbidden



#Initial Values
model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group)

# Specification data model
model.function <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  tau ~ dgamma(0.001, 0.001)
  sigma2 <- (1/tau)
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
})

samples <- nimbleMCMC(
  code = model.function,
  constants = model.data, ## provide the combined data & constants as constants
  inits = model.inits,
  monitors = parameters,
  nchains = 3,
  niter = 5000,
  nburnin = 2000,
  thin = 5)

