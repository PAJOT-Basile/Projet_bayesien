# Import libraries
require("anyLib")
anyLib(c("tidyverse", "R2jags", "lme4", "cowplot"))

# Set random seed
set.seed(2023)

# Data importation
gopher <- read.csv("gopher.csv", header=TRUE, stringsAsFactors=TRUE, sep=";", dec=".") %>% 
  mutate(standprev = (prev - mean(prev)) / sd(prev),
         Cov.y.1 = ifelse(year == 2005, 1, 0),
         Cov.y.2 = ifelse(year == 2006, 1, 0))

# The null model function
year_prevalence_model <- function(){
  # The null model does not take into account any of the variables. We just want to see the evolution of the number of shells
  # Likelihood
  for (i in 1:N){  # Loop over years
    S[i] ~ dpois(lambda[i])
    log(lambda[i]) <- mu.0 + alpha.y.1 * Cov.y.1[i] + alpha.y.2 * Cov.y.2[i] + b.prev * prev[i] + log(A[i])
  }
  # Priors
  mu.0 ~ dnorm(0, 1/100)
  alpha.y.1 ~ dnorm(0, 1/100)
  alpha.y.2 ~ dnorm(0, 1/100)
  b.prev ~ dnorm(0, 1/100)
}

# Now, we can make list of the data to use in the jags function
datax <- list(
  N = gopher$year  %>% 
    length(),
  S = gopher$shells,        # Number of shells
  A = gopher$Area,          # Area offset
  Cov.y.1 = gopher$Cov.y.1,  # Effect of the year 2005
  Cov.y.2 = gopher$Cov.y.2,  # Effect of the year 2006
  prev = gopher$prev         # Prevalence
)

# Parameters to estimate
params <- c("mu.0", "alpha.y.1", "alpha.y.2", "b.prev")

# Initialising the chains
init1 <- list(mu.0 = -0.5, alpha.y.1 = -0.5, alpha.y.2 = -0.5, b.prev = -0.5)
init2 <- list(mu.0 = 0.5, alpha.y.1 = 0.5, alpha.y.2 = 0.5, b.prev = 0.5)
init <- list(init1, init2)

# Define the iteration parameters
nb.iterations <- 9000
nb.burnin <- 4500
# Run the model using jags
M4 <- jags(
  data = datax,
  inits = init,
  parameters.to.save = params,
  model.file = year_prevalence_model,
  n.chains = 2,
  n.iter = nb.iterations,
  n.burnin = nb.burnin,
  n.thin = 1
)

M4

traceplot(M4, mfrow=c(2, 3), ask=FALSE)
par(mfrow = c(1, 1))

res <- M4$BUGSoutput$sims.matrix %>% 
  as.data.frame()
hist(res$mu.0)
hist(res$alpha.y.1)
hist(res$alpha.y.2)
hist(res$b.prev)

# The mean number of turtle shells per individual is : 
shells <- gopher$Area  * exp(res$mu.0 + res$alpha.y.1 * gopher$Cov.y.1 + res$alpha.y.2 * gopher$Cov.y.2 + res$b.prev * gopher$prev)
hist(shells)
mean(shells)
hist(gopher$shells)
