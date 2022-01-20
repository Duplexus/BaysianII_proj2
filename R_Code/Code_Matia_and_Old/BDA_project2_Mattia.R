library(nimble)
library(lattice)
library(dplyr)

setwd("~/Desktop/BDA2")

data <- readRDS("AIdataset.RDS")
head(data, 50)
str(data)
summary(data)

data[is.na(data$sofa),]   # there are 77 rows with missing values of SOFA

data$di <- data$di + 0.001  # add a small positive value to asynchrony index

data$id <- as.numeric(data$id)  # transform into numeric


#######################################################
# EXPLORATORY ANALYSIS
#######################################################

# Analysis 1

# trellis plot

xyplot(data$sofa ~ data$day|as.factor(data$id), cex=0.5,
       par.settings=list(superpose.symbol = list(pch=c(1,3,20))),
       as.table=T, auto.key=list(points=T,columns=3), xlab="day", ylab="SOFA")

# I group the ages only for the plots

data_orderage <- data %>% arrange(age)
data_orderage$classage[data_orderage$age >= 17 & data_orderage$age <= 26] <- "17-26"
data_orderage$classage[data_orderage$age >= 33 & data_orderage$age <= 44] <- "33-44"
data_orderage$classage[data_orderage$age >= 45 & data_orderage$age <= 54] <- "45-54"
data_orderage$classage[data_orderage$age >= 55 & data_orderage$age <= 64] <- "55-64"
data_orderage$classage[data_orderage$age >= 65 & data_orderage$age <= 74] <- "65-74"
data_orderage$classage[data_orderage$age >= 75 & data_orderage$age <= 80] <- "75-80"
data_orderage$classage[data_orderage$age >= 81 & data_orderage$age <= 88] <- "81-88"
data_orderage$classage <- factor(data_orderage$classage)
data_orderage$id <- factor(data_orderage$id, levels = c(unique(data_orderage$id)))

xyplot(data_orderage$sofa ~ data_orderage$day|data_orderage$id, 
       group = data_orderage$classage, cex=0.5, 
       par.settings=list(superpose.symbol = list(pch=c(1,3,20))),
       as.table=T, auto.key=list(points=T,columns=3), xlab="day", ylab="SOFA")

# Analysis 2

xyplot(data_orderage$di ~ data_orderage$day|data_orderage$id, 
       group = data_orderage$classage, cex=0.5, 
       par.settings=list(superpose.symbol = list(pch=c(1,3,20))),
       as.table=T, auto.key=list(points=T,columns=3), xlab="day", ylab="Asynchrony index")

# Analysis 3

data_orderfail <- data
data_orderfail[data_orderfail$fail==2,]$fail = 0
data_orderfail <- data_orderfail %>% arrange(fail)
data_orderfail$id <- factor(data_orderfail$id, levels = c(unique(data_orderfail$id)))

xyplot(data_orderfail$sofa ~ data_orderfail$day|data_orderfail$id, 
       group = data_orderfail$fail, cex=0.5, 
       par.settings=list(superpose.symbol = list(pch=c(1,3,20))),
       as.table=T, auto.key=list(points=T,columns=3), xlab="day", ylab="SOFA")

xyplot(data_orderfail$di ~ data_orderfail$day|data_orderfail$id, 
       group = data_orderfail$fail, cex=0.5, 
       par.settings=list(superpose.symbol = list(pch=c(1,3,20))),
       as.table=T, auto.key=list(points=T,columns=3), xlab="day", ylab="Asynchrony index")


rm(data_orderage, data_orderfail)



############################################################
# MODELS
############################################################


##########################################
# ANALYSIS 1
##########################################
# SOFA is the response variable

# NIMBLE’s default MCMC configuration will treat automatically the missing values
# as unknowns to be sampled

model.data <- list(
  y = data$sofa
)

model.constants <- list(
  x = cbind(data$day,data$age), id = data$id, N = NROW(data),  M = 139
)

###########################################################
# MODEL 1
# linear mixed model with linear covariates
###########################################################

# MODEL SPECIFICATION 
model.code <- nimbleCode({
  # Specification data model
  for (i in 1:N)
  {
    # NIMBLE allows the parametrization with the standard deviation
    y[i] ~ dnorm(y.hat[i], sd = sigma)  
    y.hat[i] <- b0[id[i]] + beta0 + beta1*x[i,1]+beta2*x[i,2]
  }
  for(j in 1:M)
  {
    b0[j] ~ dnorm(0, sd = sigma_b0)
  }
  # Prior specification
  sigma ~ dunif(0,100)
  sigma_b0 ~ dunif(0,100)
  
  beta0 ~ dnorm(0, sd = 1000)   
  beta1 ~ dnorm(0, sd = 1000)    
  beta2 ~ dnorm(0, sd = 1000)
})

# DEFINE INITIAL VALUES
model.inits <- list(
  list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = 1, sigma_b0 = 1),
  list(beta0 = 5, beta1 = -1, beta2 = -1, sigma = 5, sigma_b0 = 0.1),
  list(beta0 = -1, beta1 = 1, beta2 = 1, sigma = 0.1, sigma_b0 = 0.01)
)


model1 <- nimbleModel(code=model.code, constants=model.constants, data=model.data)

model1.conf <- configureMCMC(model1, print = TRUE)
# model1.conf$addMonitors(c("y.hat", "b0"))
model1.mcmc <- buildMCMC(model1.conf)

C.model1 <- compileNimble(model1)
C.model1.mcmc <- compileNimble(model1.mcmc, project = model1)

mcmc.out1 <- runMCMC(C.model1.mcmc, niter=100000, nburnin=10000, thin=10, nchains=3,
                     inits = model.inits, setSeed = c(3,2,1), summary = TRUE)

mcmc.out1$summary$all.chains

par(mfrow=c(2,3))
plot(mcmc.out1$samples$chain1[ , 'beta0'], type = 'l', xlab = 'iteration',  ylab = expression(beta[0]))
lines(mcmc.out1$samples$chain2[ , 'beta0'], col=2)
lines(mcmc.out1$samples$chain3[ , 'beta0'], col=3)
plot(mcmc.out1$samples$chain1[ , 'beta1'], type = 'l', xlab = 'iteration',  ylab = expression(beta[1]))
lines(mcmc.out1$samples$chain2[ , 'beta1'], col=2)
lines(mcmc.out1$samples$chain3[ , 'beta1'], col=3)
plot(mcmc.out1$samples$chain1[ , 'beta2'], type = 'l', xlab = 'iteration',  ylab = expression(beta[2]))
lines(mcmc.out1$samples$chain2[ , 'beta2'], col=2)
lines(mcmc.out1$samples$chain3[ , 'beta2'], col=3)
plot(mcmc.out1$samples$chain1[ , 'sigma'], type = 'l', xlab = 'iteration',  ylab = expression(sigma))
lines(mcmc.out1$samples$chain2[ , 'sigma'], col=2)
lines(mcmc.out1$samples$chain3[ , 'sigma'], col=3)
plot(mcmc.out1$samples$chain1[ , 'sigma_b0'], type = 'l', xlab = 'iteration',  ylab = expression(sigma[b0]))
lines(mcmc.out1$samples$chain2[ , 'sigma_b0'], col=2)
lines(mcmc.out1$samples$chain3[ , 'sigma_b0'], col=3)

###########################################################
# MODEL 2
# linear mixed model with day, day^2 and age
###########################################################

#...






#################################
# ANALYSIS 2
#################################

# Beta-regression model with random effects
# Regress the mean bounded asynchrony index as a function of day and age

model.data <- list(
  y = data$di
)

model.constants <- list(
  x = cbind(data$day,data$age), id = data$id, N = NROW(data),  M = 139
)

# MODEL SPECIFICATION 
model.code <- nimbleCode({
  for (i in 1:N)
  {
    y[i] ~ dbeta(a[i], b[i])  # standard parametrization of the beta distribution
    # reparametrizations (mu is the mean of the AI, 
    # phi can be interpreted as a precision parameter)
    a[i] <- mu[i] * phi   
    b[i] <- (1-mu[i]) * phi
    logit(mu[i]) <- b0[id[i]] + beta0 + beta1*x[i,1]+beta2*x[i,2] # regress the mean
    # phi can either have its own prior (must be greater than 0), or can be
    # modeled on covariates as well (choose a link function to keep it greater
    # than 0, such as exponential).
  }
  for(j in 1:M)
  {
    b0[j] ~ dnorm(0, sd = sigma_b0)
  }
  # Prior specification
  # phi ~ dunif(0, 500)
  # phi ~ dcauchy(0, 5)
  phi ~ dgamma(0.001, 0.001)
  sigma_b0 ~ dunif(0,100)
  
  beta0 ~ dnorm(0, sd = 1000)   
  beta1 ~ dnorm(0, sd = 1000)    
  beta2 ~ dnorm(0, sd = 1000)
})

# DEFINE INITIAL VALUES
model.inits <- list(
  list(beta0 = 0, beta1 = 0, beta2 = 0, phi = 30, sigma_b0 = 1),
  list(beta0 = 0.5, beta1 = -0.1, beta2 = -0.1, phi = 10, sigma_b0 = 1.5),
  list(beta0 = 1, beta1 = 0.1, beta2 = 0.1, phi = 20, sigma_b0 = 0.5)
)


model3 <- nimbleModel(code=model.code, constants=model.constants, data=model.data)

model3.conf <- configureMCMC(model3, print = TRUE)
model3.mcmc <- buildMCMC(model3.conf)

C.model3 <- compileNimble(model3)
C.model3.mcmc <- compileNimble(model3.mcmc, project = model3)

mcmc.out3 <- runMCMC(C.model3.mcmc, niter=100000, nburnin=30000, thin=15, nchains=3,
                     inits = model.inits, setSeed = TRUE, summary = TRUE)

mcmc.out3$summary$all.chains

par(mfrow=c(2,3))
plot(mcmc.out3$samples$chain1[ , 'beta0'], type = 'l', xlab = 'iteration',  ylab = expression(beta[0]))
lines(mcmc.out3$samples$chain2[ , 'beta0'], col=2)
lines(mcmc.out3$samples$chain3[ , 'beta0'], col=3)
plot(mcmc.out3$samples$chain1[ , 'beta1'], type = 'l', xlab = 'iteration',  ylab = expression(beta[1]))
lines(mcmc.out3$samples$chain2[ , 'beta1'], col=2)
lines(mcmc.out3$samples$chain3[ , 'beta1'], col=3)
plot(mcmc.out3$samples$chain1[ , 'beta2'], type = 'l', xlab = 'iteration',  ylab = expression(beta[2]))
lines(mcmc.out3$samples$chain2[ , 'beta2'], col=2)
lines(mcmc.out3$samples$chain3[ , 'beta2'], col=3)
plot(mcmc.out3$samples$chain1[ , 'phi'], type = 'l', xlab = 'iteration',  ylab = expression(phi))
lines(mcmc.out3$samples$chain2[ , 'phi'], col=2)
lines(mcmc.out3$samples$chain3[ , 'phi'], col=3)
plot(mcmc.out3$samples$chain1[ , 'sigma_b0'], type = 'l', xlab = 'iteration',  ylab = expression(sigma[b0]))
lines(mcmc.out3$samples$chain2[ , 'sigma_b0'], col=2)
lines(mcmc.out3$samples$chain3[ , 'sigma_b0'], col=3)
















