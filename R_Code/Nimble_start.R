library(nimble)
library(coda)
## Data
# Code from Kery and Royle (2015)
# Choose sample sizes and prepare obs. data array y
set.seed(1)                   # So we all get same data set
M <- 100                      # Number of sites
J <- 3                        # Number of repeated abundance measurements
C <- matrix(NA, nrow = M, ncol = J) # to contain the observed data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for abundance model and compute lambda
beta0 <- 0                    # Log-scale intercept
beta1 <- 2                    # Log-scale slope for vegHt
lambda <- exp(beta0 + beta1 * vegHt) # Expected abundance

# Draw local abundance
N <- rpois(M, lambda)

# Create a covariate called wind
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2                        # Logit-scale intercept
alpha1 <- -3                        # Logit-scale slope for wind
p <- plogis(alpha0 + alpha1 * wind) # Detection probability

# Take J = 3 abundance measurements at each site
for(j in 1:J) {
  C[,j] <- rbinom(M, N, p[,j])
}

# Create factors
time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = TRUE)
hab <- c(rep("A", 33), rep("B", 33), rep("C", 34))  # assumes M = 100

# Bundle data
# NIMBLE: For full flexibility, we could separate this list
#         into constants and data lists.  For simplicity we will keep
#         it as one list to be provided as the "constants" argument.
#         See comments about how we would split it if desired.
win.data <- list(
  ## NIMBLE: C is the actual data
  C = C,
  ## NIMBLE: Covariates can be data or constants
  ##         If they are data, you could modify them after the model is built
  wind = wind,
  vegHt = vegHt,
  XvegHt = seq(-1, 1,, 100), # Used only for derived quantities
  Xwind = seq(-1, 1,,100),   # Used only for derived quantities
  ## NIMBLE: The rest of these are constants, needed for model definition
  ## We can provide them in the same list and NIMBLE will figure it out.
  M = nrow(C),
  J = ncol(C),
  hab = as.numeric(factor(hab))
)
#### Model ####
Section6p4_code <- nimbleCode( {
  # Priors
  for(k in 1:3) {                # Loop over 3 levels of hab or time factors
    alpha0[k] ~ dunif(-10, 10) # Detection intercepts
    alpha1[k] ~ dunif(-10, 10) # Detection slopes
    beta0[k] ~ dunif(-10, 10)  # Abundance intercepts
    beta1[k] ~ dunif(-10, 10)  # Abundance slopes
  }
  
  # Likelihood
  # Ecological model for true abundance
  for (i in 1:M){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0[hab[i]] + beta1[hab[i]] * vegHt[i]
    # Some intermediate derived quantities
    critical[i] <- step(2-N[i])# yields 1 whenever N is 2 or less
    z[i] <- step(N[i]-0.5)     # Indicator for occupied site
    # Observation model for replicated counts
    for (j in 1:J){
      C[i,j] ~ dbin(p[i,j], N[i])
      logit(p[i,j]) <- alpha0[j] + alpha1[j] * wind[i,j]
    }
  }
  
  # Derived quantities; unnececssary when running for inference purpose
  # NIMBLE: We have filled in indices in the next two lines.
  Nocc <- sum(z[1:100])         # Number of occupied sites among sample of M
  Ntotal <- sum(N[1:100])       # Total population size at M sites combined
  Nhab[1] <- sum(N[1:33])  # Total abundance for sites in hab A
  Nhab[2] <- sum(N[34:66]) # Total abundance for sites in hab B
  Nhab[3] <- sum(N[67:100])# Total abundance for sites in hab C
  for(k in 1:100){         # Predictions of lambda and p ...
    for(level in 1:3){    #    ... for each level of hab and time factors
      lam.pred[k, level] <- exp(beta0[level] + beta1[level] * XvegHt[k])
      logit(p.pred[k, level]) <- alpha0[level] + alpha1[level] * Xwind[k]
    }
  }
  # NIMBLE: We have filled in indices in the next line. 
  N.critical <- sum(critical[1:100]) # Number of populations with critical size
})

Nst <- apply(C, 1, max)+1   # Important to give good inits for latent N
inits <- function() list(N = Nst, 
                         alpha0 = rnorm(3), 
                         alpha1 = rnorm(3), 
                         beta0 = rnorm(3), 
                         beta1 = rnorm(3))

# Parameters monitored
# could also estimate N, bayesian counterpart to BUPs before: simply add "N" to the list
params <- c("alpha0", "alpha1", "beta0", "beta1", "Nocc", "Ntotal", "Nhab", "N.critical", "lam.pred", "p.pred")


samples <- nimbleMCMC(
  code = Section6p4_code,
  constants = win.data, ## provide the combined data & constants as constants
  inits = inits,
  monitors = params,
  niter = 3000,
  nburnin = 2000,
  thin = 1)
coda.samples <- as.mcmc(samples)