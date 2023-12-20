### Model with discontinuous prevalence ### 

## load libraries ##
library(tidyverse)
library(cowplot)
library(R2jags)

## import dataset
gopher <- read.csv(file = "./gopher.csv", header = TRUE, sep = ";", dec=".", stringsAsFactors = TRUE)
gopher$year<-as.factor(gopher$year) # set year as factor
gopher$standprev<-(gopher$prev-mean(gopher$prev))/sd(gopher$prev) #standardize prev 
gopher$H<-ifelse(gopher$prev>= 25,1,0) #tranform prev as discontinuous with H=0 if prev<25% and H=1 if prev>=25%


set.seed(2023)
gopher_discont_prev <-function(){
  for (i in 1:N){
    S[i] ~ dpois(lambda[i])
    log(lambda[i]) <- mu.0 + b.prev * disc_prev[i] + log(A[i])
  }
  mu.0 ~ dnorm(0,0.001)
  b.prev ~ dnorm(0,0.001)
}

datax <- list(N = length(gopher$shells),
              disc_prev = gopher$H,
              S = gopher$shells,
              A = gopher$Area
)
init1 <- list(mu.0=0.5,b.prev=0.5)
init2 <- list(mu.0= -0.5,b.prev= -0.5)
inits <- list(init1,init2)

params <- c("mu.0","b.prev")

M7 <- jags(data=datax,
           inits=inits,
           parameters.to.save = params,
           model.file=gopher_discont_prev,
           n.chains = 2,
           n.iter = 9000,
           n.burnin = 4500,
           n.thin=1)

M7

#
traceplot(M7,mfrow=c(3,1),varname=c('mu.0','b.prev','deviance'),ask=FALSE)
head(M7$BUGSoutput$sims.matrix)
acf(M7$BUGSoutput$sims.matrix[,2])
