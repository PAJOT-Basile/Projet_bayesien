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
  for (i in 1:n){
    y[i]~dpois(lambda[i])
    log(lambda[i])<-mu + beta*x[i] + log(A[i])
  }
  mu~dnorm(0,0.001)
  beta~dnorm(0,0.001)
}

datax <- list(n=length(gopher$shells),
              x=gopher$H,
              y=gopher$shells,
              A = gopher$Area
)
init1 <- list(mu=0.5,beta=0.5)
init2 <- list(mu= -0.5,beta= -0.5)
inits <- list(init1,init2)

params <- c("mu","beta")

M2 <- jags(data=datax,
           inits=inits,
           parameters.to.save = params,
           model.file=gopher_discont_prev,
           n.chains = 2,
           n.iter = 10000,
           n.burnin = 100,
           n.thin=1)

M2

#
traceplot(M2,mfrow=c(3,1),varname=c('mu','beta','deviance'),ask=FALSE)
head(M2$BUGSoutput$sims.matrix)
acf(M2$BUGSoutput$sims.matrix[,2])
