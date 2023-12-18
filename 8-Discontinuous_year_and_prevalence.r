# Import libraries
require("anyLib")
anyLib(c("tidyverse", "R2jags", "lme4", "cowplot"))

# Set random seed
set.seed(2023)

# Data importation
gopher <- read.csv("gopher.csv", header=TRUE, stringsAsFactors=TRUE, sep=";", dec=".") %>% 
  mutate(H = ifelse(prev >= 25, 1, 0),
         Cov.y.1 = ifelse(year == 2005, 1, 0),
         Cov.y.2 = ifelse(year == 2006, 1, 0))

# The random model
disc_year_prev_model <- function(){
  # This model takes into account a random effect for the site
  # Likelihood
  for(i in 1:N){
    S[i] ~ dpois(lambda[i])
    log(lambda[i]) <- mu.0 + b.prev * prev[i] + alpha.y.1 * Cov.y.1[i] + alpha.y.2 * Cov.y.2[i] + log(A[i])
  }
  # Priors
  mu.0 ~ dnorm(0, 1/100)
  b.prev ~ dnorm(0, 1/100)
  alpha.y.1 ~ dnorm(0, 1/100)
  alpha.y.2 ~ dnorm(0, 1/100)
}

# Make the data to use in jags
datax <- list(
  N = gopher$year %>% 
    length(),
  S = gopher$shells,
  prev = gopher$H,
  A = gopher$Area,
  Cov.y.1 = ifelse(gopher$year == 2005, 1, 0),
  Cov.y.2 = ifelse(gopher$year == 2006, 1, 0)
)

# Make a list of parameters to save
params = c("mu.0", "b.prev", "alpha.y.1", "alpha.y.2")

# Initial conditions
init1 <- list(
  "mu.0" = 0.5,
  "b.prev" = 0.5,
  "alpha.y.1" = 0.5,
  "alpha.y.2" = 0.5
)
init2 <- list(
  "mu.0" = - 0.5,
  "b.prev" = - 0.5,
  "alpha.y.1" = -0.5,
  "alpha.y.2" = -0.5
)
init <- list(init1, init2)

# Iteration parameters
nb.iterations <- 9000
nb.burnin <- 4500

# Run the model
M8 <- jags(
  data = datax,
  parameters.to.save = params,
  inits = init,
  model.file = disc_year_prev_model,
  n.chains = 2,
  n.iter = nb.iterations,
  n.burnin = nb.burnin,
  n.thin = 1
)

M8
traceplot(M8, mfrow=c(2, 3), ask=FALSE)
par(mfrow=c(1, 1))

res <- M8$BUGSoutput$sims.matrix %>% 
  as.data.frame()

# Looking at the distribution
hist(res$b.prev)
hist(res$mu.0)
hist(res$alpha.y.1)
hist(res$alpha.y.2)

# Calculating the mean number of shells per individual 
shells <- gopher$Area * exp(res$mu.0 + res$b.prev * gopher$H + res$alpha.y.1 * gopher$Cov.y.1 + res$alpha.y.2 * gopher$Cov.y.2)

hist(shells)
mean(shells)
