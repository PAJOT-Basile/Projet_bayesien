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
null_year_site_model <- function(){
  # The null model does not take into account any of the variables. We just want to see the evolution of the number of shells
  # Likelihood
  for (i in 1:N){  # Loop over years
    S[i] ~ dpois(lambda[i])
    log(lambda[i]) <- mu.0 + alpha.y.1 * Cov.y.1[i] + alpha.y.2 * Cov.y.2[i] + alpha.s[site[i]] + log(A[i])
  }
  # Priors
  mu.0 ~ dnorm(0, 1/100)
  alpha.y.1 ~ dnorm(0, 1/100)
  alpha.y.2 ~ dnorm(0, 1/100)
  for (j in 1:nb.sites){
    alpha.s[j] ~ dnorm(0, tau.s)
  }
  tau.s <- 1/(sigma.s * sigma.s)
  sigma.s ~ dunif(0, 100)
}

# Now, we can make list of the data to use in the jags function
datax <- list(
  N = gopher$year  %>% 
    length(),
  S = gopher$shells,        # Number of shells
  A = gopher$Area,          # Area offset
  Cov.y.1 = gopher$Cov.y.1,  # Effect of the year 2005
  Cov.y.2 = gopher$Cov.y.2,   # Effect of the year 2006
  nb.sites = gopher$Site %>%  # Number of sites
    unique() %>% 
    length(),
  site = gopher$Site %>%    # Vector of different sites as different numbers
    as.numeric()
)

# Parameters to estimate
params <- c("mu.0", "alpha.y.1", "alpha.y.2", "sigma.s")

# Initialising the chains
init1 <- list(mu.0 = -0.5, alpha.y.1=-0.5, alpha.y.2=-0.5, sigma.s = 0.5)
init2 <- list(mu.0 = 0.5, alpha.y.1=0.5, alpha.y.2=0.5, sigma.s = 0.25)
init <- list(init1, init2)

# Define the iteration parameters
nb.iterations <- 9000
nb.burnin <- 4500
# Run the model using jags
M5 <- jags(
  data = datax,
  inits = init,
  parameters.to.save = params,
  model.file = null_year_site_model,
  n.chains = 2,
  n.iter = nb.iterations,
  n.burnin = nb.burnin,
  n.thin = 1
)

M5

traceplot(M5, mfrow=c(2, 3), ask=FALSE)
par(mfrow = c(1, 1))

res <- M5$BUGSoutput$sims.matrix %>% 
  as.data.frame()
hist(res$mu.0)
hist(res$alpha.y.1)
hist(res$alpha.y.2)
hist(res$sigma.s)

# The mean number of turtle shells per individual is : 
shells <- gopher$Area  * exp(res$mu.0 + res$alpha.y.1 * gopher$Cov.y.1 + res$alpha.y.2 * gopher$Cov.y.2 + rnorm(1, mean=0, sd=res$sigma.s))
hist(shells)
mean(shells)
hist(gopher$shells)
