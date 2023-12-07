# Import libraries
require("anyLib")
anyLib(c("tidyverse", "R2jags", "lme4", "cowplot"))

# Set random seed
set.seed(2023)

# Data importation
gopher <- read.csv("gopher.csv", header=TRUE, stringsAsFactors=TRUE, sep=";", dec=".") %>% 
  mutate(standprev = (prev - mean(prev)) / sd(prev))

# The random model
random_model <- function(){
  # This model takes into account a random effect for the site
  # Likelihood
  for(i in 1:N){
    S[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha.s[site[i]] + b.prev * prev[i] + mu.0 + alpha.y[year[i]] + log(A[i])
  }
  # Priors
  mu.0 ~ dnorm(0, 1/100)
  b.prev ~ dnorm(0, 1/100)
  for (k in 1:nb.years){
    alpha.y[k] ~ dnorm(0, 1/100)
  }
  for (j in 1:nb.sites){
    alpha.s[j] ~ dnorm(mu.s, tau.s)
  }
  mu.s ~ dnorm(0, 1/100)
  tau.s <- 1 / (sd.s * sd.s)
  sd.s ~ dunif(0, 100)
}

# Make the data to use in jags
datax <- list(
  N = gopher$year %>% 
    length(),
  S = gopher$shells,
  prev = gopher$standprev,
  A = gopher$Area,
  site = gopher$Site %>% 
    as.numeric(),
  nb.sites = gopher$Site %>% 
    unique() %>% 
    length(),
  year = ifelse(gopher$year == 2004, 1, ifelse(gopher$year == 2005, 2, 3)),
  nb.years = gopher$year %>% 
    unique() %>% 
    length()
)

# Make a list of parameters to save
params = c("mu.0", "b.prev", "mu.s", "sd.s", "alpha.y")

# Initial conditions
init1 <- list(
  "mu.0" = 0.5,
  "b.prev" = 0.5,
  "mu.s" = 0.5,
  "sd.s" = 0.5,
  "alpha.y" = rep(0.5, gopher$year %>% unique() %>% length())
)
init2 <- list(
  "mu.0" = - 0.5,
  "b.prev" = - 0.5,
  "mu.s" = - 0.5,
  "sd.s" =  0.25,
  "alpha.y" = rep(-0.5, gopher$year %>% unique() %>% length())
)
init <- list(init1, init2)

# Iteration parameters
nb.iterations <- 100000
nb.burnin <- 1000

# Run the model
M2 <- jags(
  data = datax,
  parameters.to.save = params,
  inits = init,
  model.file = random_model,
  n.chains = 2,
  n.iter = nb.iterations,
  n.burnin = nb.burnin,
  n.thin = 1
)

M2
traceplot(M2, mfrow=c(2, 3), ask=FALSE)
par(mfrow=c(1, 1))

autocorr.plot(as.mcmc(M2), ask=FALSE)
par(mfrow=c(1, 1))
