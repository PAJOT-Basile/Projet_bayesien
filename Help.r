# Import libraries
require("anyLib")
anyLib(c("tidyverse", "R2jags", "lme4", "cowplot"))

# Set random seed
set.seed(2023)

# Data importation
gopher <- read.csv("gopher.csv", header=TRUE, stringsAsFactors=TRUE, sep=";", dec=".") %>% 
  mutate(standprev = (prev - mean(prev)) / sd(prev),
         H = ifelse(prev <= 25, 0, 1),
         Cov.y.1 = ifelse(year == 2005, 1, 0),
         Cov.y.2 = ifelse(year == 2006, 1, 0))

random_model <- function(){
  # This model takes into account a random effect for the site
  # Likelihood
  for(i in 1:N){
    S[i] ~ dpois(lambda[i])
    log(lambda[i]) <- mu.0 + alpha.s[site[i]] + b.prev * prev[i] + alpha.y.1 * Cov.y.1[i] + alpha.y.2 * Cov.y.2[i] + log(A[i])
  }
  # Priors
  mu.0 ~ dnorm(0, 1/100)
  b.prev ~ dnorm(0, 1/100)
  for (j in 1:nb.sites){
    alpha.s[j] ~ dnorm(0, tau.s)
  }
  tau.s <- 1 / (sd.s * sd.s)
  sd.s ~ dunif(0, 100)
  alpha.y.1 ~ dnorm(0, 1/100)
  alpha.y.2 ~ dnorm(0, 1/100)
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
  Cov.y.1 = ifelse(gopher$year == 2005, 1, 0),
  Cov.y.2 = ifelse(gopher$year == 2006, 1, 0)
)

# Make a list of parameters to save
params = c("mu.0", "b.prev", "sd.s", "alpha.y.1", "alpha.y.2")

# Initial conditions
init1 <- list(
  "mu.0" = 0.5,
  "b.prev" = 0.5,
  "sd.s" = 0.5,
  "alpha.y.1" = 0.5,
  "alpha.y.2" = 0.5
)
init2 <- list(
  "mu.0" = - 0.5,
  "b.prev" = - 0.5,
  "sd.s" =  0.25,
  "alpha.y.1" = -0.5,
  "alpha.y.2" = -0.5
)
init <- list(init1, init2)

# Iteration parameters
nb.iterations <- 9000
nb.burnin <- 4500

# Run the model
M6 <- jags(
  data = datax,
  parameters.to.save = params,
  inits = init,
  model.file = random_model,
  n.chains = 2,
  n.iter = nb.iterations,
  n.burnin = nb.burnin,
  n.thin = 1
)

# Regardons le modèle et les traces
M6
traceplot(M6, mfrow=c(2, 3), ask=FALSE)
par(mfrow=c(1, 1))

# On récupère les paramètres et on regarde leur distribution
res <- M6$BUGSoutput$sims.matrix %>% 
  as.data.frame()
hist(res$b.prev)
hist(res$mu.0)
hist(res$alpha.y.1)
hist(res$alpha.y.2)
hist(res$sd.s)

# Rétrotransformation : idée 1
# Moyennes
beta <- mean(res$b.prev)
alpha.1 <- mean(res$alpha.y.1)
alpha.2 <- mean(res$alpha.y.2)
mu.0 <- mean(res$mu.0)


shells6 <- matrix(NA, ncol=nrow(gopher), nrow=nrow(res))
for (i in 1:nrow(gopher)){
  shells6[, i] <- gopher$Area[i] * exp(res$mu.0 + 
                                         res$b.prev * gopher$prev[i] + 
                                         res$alpha.y.1 * gopher$Cov.y.1[i] + 
                                         res$alpha.y.2 * gopher$Cov.y.2[i] + 
                                         rnorm(nrow(res), mean=0, sd=res$sd.s)) 
}
head(shells6)
hist(shells6)
mean(shells6)








simulated_shells <- rep(NA, nrow(res))
for (i in 1:nrow(res)){
  lambda <- gopher$Area * exp(res$mu.0[i] + res$alpha.y.1[i] * gopher$Cov.y.1 + res$alpha.y.2[i] * gopher$Cov.y.2 + res$b.prev[i] * gopher$standprev + rnorm(1, mean=0, sd=res$sd.s[i]))
  simulated_shells[i] <- rpois(1, lambda)
}
hist(simulated_shells)
mean(simulated_shells)

# Rétrotransformation : idée 2
simulated_shells <- gopher$Area * exp(res$mu.0 + res$alpha.y.1 * gopher$Cov.y.1 + res$alpha.y.2 * gopher$Cov.y.2 + res$b.prev * gopher$standprev + rnorm(1, mean=0, sd=res$sd.s)) 
hist(simulated_shells)
mean(simulated_shells)

# Récupérons le DIC
DIC6 <- M6$BUGSoutput$DIC