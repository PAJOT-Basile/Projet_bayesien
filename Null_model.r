# Import libraries
require("anyLib")
anyLib(c("tidyverse", "R2jags", "lme4", "cowplot"))

# Set random seed
set.seed(2023)

# Data importation
gopher <- read.csv("gopher.csv", header=TRUE, stringsAsFactors=TRUE, sep=";", dec=".") %>% 
  mutate(standprev = (prev - mean(prev)) / sd(prev))

# The null model function
null_model <- function(){
  # The null model does not take into account any of the variables. We just want to see the evolution of the number of shells
  # Likelihood
  for (i in 1:N){  # Loop over years
    S[i] ~ dpois(lambda[i])
    log(lambda[i]) <- mu.0 + log(A[i])
  }
  # Priors
  mu.0 ~ dnorm(0, 1/100)
}

# Now, we can make list of the data to use in the jags function
datax <- list(
  N = gopher$year  %>% 
    length(),
  S = gopher$shells,        # Number of shells
  A = gopher$Area          # Area offset
  #d = gopher$density        # density of turtoise
)

# Parameters to estimate
params <- c("mu.0")

# Initialising the chains
init1 <- list(mu.0 = -0.5)
init2 <- list(mu.0 = 0.5)
init <- list(init1, init2)

# Define the iteration parameters
nb.iterations <- 10000
nb.burnin <- 100
# Run the model using jags
M0 <- jags(
  data = datax,
  inits = init,
  parameters.to.save = params,
  model.file = null_model,
  n.chains = 2,
  n.iter = nb.iterations,
  n.burnin = nb.burnin,
  n.thin = 1
)

M0

traceplot(M0, mfrow=c(1, 2), ask=FALSE)
par(mfrow = c(1, 1))

res <- M0$BUGSoutput$sims.matrix %>% 
  as.data.frame()
hist(res$mu.0)

# The mean number of turtle shells per individual is : 
shells <- gopher$Area  * exp(res$mu.0)
hist(shells)
mean(shells)
hist(gopher$shells)
