library(tidyverse)
library(sf)
library(sp)
library(raster)
library(rjags)
library(runjags)
library(coda)
library(ggmcmc)
library(MCMCvis)
library(arm)
library(loo)


#Bring in data
model.data <- read.csv('model_data.csv')


#Center/scale covariates
mean.prop <- mean(model.data$prop)
prop2sd <- 2*sd(model.data$prop)
model.data$cs.prop <- NA
for(i in 1:length(model.data$prop)){
  model.data$cs.prop[i] <- (model.data$prop[i] - mean.prop)/prop2sd
}

model.data$prop.sq <- rescale(I(model.data$prop^2))

n.db <- length(unique(model.data$DetNum))
n.samps <- length(model.data$redetects)


##JAGS code##

#psi is the proportion of redetects > 0
#model code below allows for switching between Poisson regression and negative binomial
#data not overdispersed
#commented section to account for calving grounds since behavior and effort are different
#Iterate over each resurvey, i

##Intercept-only model##
effort_mod.int <- "model{
#Likelihood statement
for(i in 1:n.samps){

#Account for zero inflation
mu[i] <- lambdaN[i]*z[i] + 0.00001 #uncomment for Poisson regression
#mu[i] <- lambda.star[i]*z[i] + 0.00001 #uncomment for negative binomial
z[i] ~ dbern(psi)

N[i] ~ dpois(mu[i])
#lambda.star[i] <- lambdaN[i] * r[i]
#r[i] ~ dgamma(r.N, r.N)

lambdaN[i] <- exp(alpha)

loglik[i] <- log(dpois(N[i],mu[i])) 

#Create replicate redetections (for Bayesian p-values)
Nnew[i]~dpois(mu[i]) 

#Create fit statistic 1
	FT1[i]<-pow(sqrt(N[i])-sqrt(lambdaN[i]),2)
	FT1new[i]<-pow(sqrt(Nnew[i])-sqrt(lambdaN[i]),2)

}

T1p <- sum(FT1[1:n.samps])
T1newp <- sum(FT1new[1:n.samps])
Bp.N <- sum(T1newp) > sum(T1p)


#Priors
alpha ~ dnorm(0, 0.1)
psi ~ dunif(0, 1)
#r.N ~ dunif(0, 100)

}"


params <- c("alpha", "psi", "Bp.N", "loglik") #"r.N",

inits <- function(){ list(alpha = rnorm(1, 0, 0.1),
                          psi = runif(1, 0, 1))
}

jags.dat <- list(n.samps = n.samps,
                 N = model.data$redetects)

rd.per.e.int <- run.jags(model = effort_mod.int,
                     monitor = params,
                     data = jags.dat,
                     n.chains = 3,
                     burnin = 5000,
                     sample = 8000,
                     inits = inits) 

#Diagnostics
#Easier to just run the model with/without loglik than to remove that parameter
effort.mod_mcmc <- as.mcmc.list(rd.per.e.int)
effort.mod_ggs <- ggs(effort.mod_mcmc)
effort.mod_ggs <- effort.mod_ggs %>% filter(!grepl("loglik", Parameter))
ggs_geweke(effort.mod_ggs)
ggs_Rhat(effort.mod_ggs)

ggs_traceplot(effort.mod_ggs, c("alpha"))
ggs_traceplot(effort.mod_ggs, c("psi"))

effort.mod.summ.int <- MCMCsummary(effort.mod_mcmc)


##Quadratic model##
effort_mod.q <- "model{
#Likelihood statement
for(i in 1:n.samps){

#Account for zero inflation
mu[i] <- lambdaN[i]*z[i] + 0.00001 #uncomment for Poisson regression
#mu[i] <- lambda.star[i]*z[i] + 0.00001 #uncomment for negative binomial
z[i] ~ dbern(psi)

N[i] ~ dpois(mu[i])
#lambda.star[i] <- lambdaN[i] * r[i]
#r[i] ~ dgamma(r.N, r.N)

##Use options below to include/exclude calf parameter
#lambdaN[i] <- exp(alpha + b.prop * prop[i] + b.prop.sq * prop.sq[i] + eps.db[DetNum[i]])
lambdaN[i] <- exp(alpha + b.prop * prop[i] + b.prop.sq * prop.sq[i] + b.calf*calf[i] + eps.db[DetNum[i]])

loglik[i] <- log(dpois(N[i],mu[i])) 

#Create replicate redetections (for Bayesian p-values)
Nnew[i]~dpois(mu[i]) 

#Create fit statistic 1
	FT1[i]<-pow(sqrt(N[i])-sqrt(lambdaN[i]),2)
	FT1new[i]<-pow(sqrt(Nnew[i])-sqrt(lambdaN[i]),2)

}

T1p <- sum(FT1[1:n.samps])
T1newp <- sum(FT1new[1:n.samps])
Bp.N <- sum(T1newp) > sum(T1p)

#Account for DetNum as random effect
for(j in 1:max(DetNum)){
eps.db[j] ~ dnorm(0, tau.eps.db)
}

#Priors
alpha ~ dnorm(0, 0.1)
b.prop ~ dnorm(0, 0.1)
b.prop.sq ~ dnorm(0, 0.1)
b.calf ~ dnorm(0, 0.1) #comment out for non-calf model
psi ~ dunif(0, 1)
tau.eps.db ~ dgamma(1, 1)
sigma.eps.db <- 1/sqrt(tau.eps.db)
#r.N ~ dunif(0, 100)

}"


params <- c("alpha", "b.prop", "b.prop.sq", "sigma.eps.db", "psi", "Bp.N", "loglik", "b.calf") #"r.N"

inits <- function(){ list(alpha = rnorm(1, 0, 0.1),
                          b.prop = rnorm(1, 0, 0.1),
                          b.prop.sq = rnorm(1, 0, 0.1),
                          b.calf = rnorm(1, 0, 0.1), #comment out for non-calf model
                          psi = runif(1, 0, 1))
}

jags.dat <- list(n.samps = n.samps,
                 prop = model.data$cs.prop,
                 prop.sq = model.data$prop.sq,
                 calf = model.data$calf, #comment out for non-calf model
                 N = model.data$redetects,
                 DetNum = model.data$DetNum)

rd.per.e.sq <- run.jags(model = effort_mod.q,
                         monitor = params,
                         data = jags.dat,
                         n.chains = 3,
                         burnin = 5000,
                         sample = 8000,
                         inits = inits) 

#Diagnostics
effort.mod_mcmc <- as.mcmc.list(rd.per.e.sq)
effort.mod_ggs <- ggs(effort.mod_mcmc)
effort.mod_ggs <- effort.mod_ggs %>% filter(!grepl("loglik", Parameter))
ggs_geweke(effort.mod_ggs)
ggs_Rhat(effort.mod_ggs)

ggs_traceplot(effort.mod_ggs, c("alpha"))
ggs_traceplot(effort.mod_ggs, c("b"))
ggs_traceplot(effort.mod_ggs, c("psi"))

effort.mod.summ.sq <- MCMCsummary(effort.mod_mcmc)



##Linear model##
effort_mod.l <- "model{
#Likelihood statement
for(i in 1:n.samps){

#Account for zero inflation
mu[i] <- lambdaN[i]*z[i] + 0.00001 #uncomment for Poisson regression
#mu[i] <- lambda.star[i]*z[i] + 0.00001 #uncomment for negative binomial
z[i] ~ dbern(psi)

N[i] ~ dpois(mu[i])
#lambda.star[i] <- lambdaN[i] * r[i]
#r[i] ~ dgamma(r.N, r.N)

#lambdaN[i] <- exp(alpha + b.prop * prop[i] + eps.db[DetNum[i]])
lambdaN[i] <- exp(alpha + b.prop * prop[i] + b.calf*calf[i]  + eps.db[DetNum[i]])

loglik[i] <- log(dpois(N[i],mu[i]))

#Create replicate redetections (for Bayesian p-values)
Nnew[i]~dpois(mu[i])

#Create fit statistic 1
	FT1[i]<-pow(sqrt(N[i])-sqrt(lambdaN[i]),2)
	FT1new[i]<-pow(sqrt(Nnew[i])-sqrt(lambdaN[i]),2)

}

T1p <- sum(FT1[1:n.samps])
T1newp <- sum(FT1new[1:n.samps])
Bp.N <- sum(T1newp) > sum(T1p)

#Account for DetNum as random effect
for(j in 1:max(DetNum)){
eps.db[j] ~ dnorm(0, tau.eps.db)
}

#Priors
alpha ~ dnorm(0, 0.1)
b.prop ~ dnorm(0, 0.1)
b.calf ~ dnorm(0, 0.1) #comment out for non-calf model
psi ~ dunif(0, 1)
tau.eps.db ~ dgamma(1, 1)
sigma.eps.db <- 1/sqrt(tau.eps.db)
#r.N ~ dunif(0, 100)

}"


params <- c("alpha", "b.prop", "sigma.eps.db", "psi", "Bp.N", "loglik", "b.calf") #"r.N", 

inits <- function(){ list(alpha = rnorm(1, 0, 0.1),
                          b.prop = rnorm(1, 0, 0.1),
                          b.calf = rnorm(1, 0, 0.1), #comment out for non-calf model
                          psi = runif(1, 0, 1))
}

jags.dat <- list(n.samps = n.samps,
                 prop = model.data$cs.prop,
                 calf = model.data$calf, #comment out for non-calf model
                 N = model.data$redetects,
                 DetNum = model.data$DetNum)

rd.per.e.l <- run.jags(model = effort_mod.l,
                     monitor = params,
                     data = jags.dat,
                     n.chains = 3,
                     burnin = 5000,
                     sample = 8000,
                     inits = inits) 

#Run diagnostics 
effort.mod_mcmc <- as.mcmc.list(rd.per.e.l)
effort.mod_ggs <- ggs(effort.mod_mcmc)
effort.mod_ggs <- effort.mod_ggs %>% filter(!grepl("loglik", Parameter))
ggs_geweke(effort.mod_ggs)
ggs_Rhat(effort.mod_ggs)

ggs_traceplot(effort.mod_ggs, c("alpha"))
ggs_traceplot(effort.mod_ggs, c("b"))
ggs_traceplot(effort.mod_ggs, c("psi"))

effort.mod.summ.l <- MCMCsummary(effort.mod_mcmc)


#Compare models
intercept.results <- as.mcmc.list(rd.per.e.int)
linear.results <- as.mcmc.list(rd.per.e.l)
quadratic.results <- as.mcmc.list(rd.per.e.sq)
loglikint <- combine.mcmc(intercept.results, vars = c("loglik"))
loglikl <- combine.mcmc(linear.results, vars = c("loglik"))
logliksq <- combine.mcmc(quadratic.results, vars = c("loglik"))
# store the log likelihood values as a matrix
logliki.mat <- as.matrix(loglikint)
loglikl.mat <- as.matrix(loglikl)
logliksq.mat <- as.matrix(logliksq)

waic(logliki.mat)
loo(logliki.mat)
waic(loglikl.mat)
loo(loglikl.mat)
waic(logliksq.mat)
loo(logliksq.mat)


