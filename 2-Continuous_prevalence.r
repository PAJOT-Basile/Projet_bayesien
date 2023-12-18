# Import libraries
require("anyLib")
anyLib(c("tidyverse", "R2jags", "lme4", "cowplot"))

# Set random seed
set.seed(2023)

# Data importation
gopher <- read.csv("gopher.csv", header=TRUE, stringsAsFactors=TRUE, sep=";", dec=".") %>% 
  mutate(standprev = (prev - mean(prev)) / sd(prev))

# The prevalence model
prevalence_model <- function(){
  # We here use the null model and add the prevalence and a its coefficient
  # Likelihood
  for (i in 1:N){
    S[i] ~ dpois(lambda[i])
    log(lambda[i]) <- mu.0 + b.prev * prev[i] + log(A[i])
  }
  # Priors
  mu.0 ~ dnorm(0, 1/100)
  b.prev ~ dnorm(0, 1/100)
}

# We make a list of the data to use in jags
datax <- list(
  N = gopher$year %>% 
    length(),
  S = gopher$shells,
  prev = gopher$standprev,
  A = gopher$Area
)

# Parameters to estimate
params <- c("mu.0", "b.prev")

# We initialise the model
init1 <- list("mu.0"=0.5, "b.prev"=0.5)
init2 <- list("mu.0"=-0.5, "b.prev"=-0.5)
init <- list(init1, init2)

# We define the iteration parameters
nb.iterations <- 9000
nb.burnin <- 4500

# We run the model
M2 <- jags(
  data = datax,
  inits = init,
  parameters.to.save = params,
  model.file = prevalence_model,
  n.chains = 2,
  n.iter = nb.iterations,
  n.burnin = nb.burnin,
  n.thin=1
)

M2

traceplot(M2, mfrow=c(1, 3), ask=FALSE)
par(mfrow=c(1, 1))

res <- M2$BUGSoutput$sims.matrix %>% 
  as.data.frame()

# Looking at the distribution
hist(res$b.prev)
hist(res$mu.0)

# Calculating the mean number of shells per individual 
shells <- gopher$Area * exp(res$mu.0 + res$b.prev*gopher$standprev)

hist(shells)
mean(shells)
