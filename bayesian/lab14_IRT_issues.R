################################################################
## Bayesian Statistics: Advanced IRT
## Quant III Lab 13
## December 5th 2013
################################################################

# install.packages("wnominate")
# install.packages("msm")
# install.packages("pscl")

library(pscl)
library(rstan)

# loading roll call data
# Source: http://jackman.stanford.edu/blog/
load("lab13_senate_rc.rda")
rc <- dropRollCall(rc, dropList=list(codes = "notInLegis", lop = 0))


################################################################
## BASELINE IRT MODEL
################################################################

irt <- ideal(rc, store.item=TRUE)

# analysis of results
summary(irt)

################################################################
## ASSESSING CONVERGENCE
################################################################

# Function to convert ideal objects to coda/mcmc objects
legislators_to_coda <- function(irt) mcmc(irt$x[,,1])
legislators_to_array <- function(irt) aperm(irt$x, c(1,3,2))
items_to_coda <- function(irt, par) mcmc(irt$beta[,,par])
items_to_array <- function(irt, par){
    t <- irt$beta[,,par]
    a <- array(t, dim=c(dim(t)[1], 1, dim(t)[2]),
        dimnames=list(dimnames(t)[[1]], 1, dimnames(t)[[2]]))
    return(a)
}


# Visual test
plot(legislators_to_coda(irt), ask=TRUE)
plot(items_to_coda(irt, 'Discrimination D1'), ask=TRUE)

# Summary with rstan
monitor(legislators_to_array(irt))
monitor(items_to_array(irt, 'Discrimination D1'))

### SOLUTIONS (?)

# If chain has not converged, longer chain can work
irt <- ideal(rc, store.item=TRUE, maxiter=50000, thin=200, burnin=10000,
    verbose=TRUE)

load("lab14_irt_long_chain.Rdata")

# Convergence
plot(legislators_to_coda(irt), ask=TRUE)
plot(items_to_coda(irt, 'Discrimination D1'), ask=TRUE)
monitor(legislators_to_array(irt))
monitor(items_to_array(irt, 'Discrimination D1'))

# Another alternative is to use a hierarchical approach
irt <- ideal(rc, store.item=TRUE, normalize=TRUE)
# 'normalize' identifies the model imposing the constraint that ideal points
# have unit variance in each dimension (mean 0 and sd 1)
# This is equivalent to a hierarchical model where x_i ~ N(0, 1)
# where x_i has an informative prior distribution (we fix the
# hyperparameters)

monitor(legislators_to_array(irt))
monitor(items_to_array(irt, 'Discrimination D1'))

# Other solution: parameter expansion (see BDA for more details)
irt <- ideal(rc, store.item=TRUE, mda=TRUE)

################################################################
## IRT WITH >1 DIMENSIONS
################################################################

# 'ideal' function works with multiple dimensions
irt <- ideal(rc, d=2, store.item=TRUE, maxiter=50000, thin=200, 
    burnin=10000, verbose=TRUE)

load("lab14_irt_2D.Rdata")

# low values in first dimension
head(irt$xbar[order(irt$xbar[,1]),])
# high values in first dimension
tail(irt$xbar[order(irt$xbar[,1]),])

# low values in second dimension
head(irt$xbar[order(irt$xbar[,2]),])
# high values in second dimension
tail(irt$xbar[order(irt$xbar[,2]),])

# items that discriminate in second dimension
discrimination <- irt$betabar[,"Discrimination D2"]

# top 2 most "discriminatory" bills for POSITIVE values of scale
rc$vote.data[order(discrimination, decreasing=TRUE)[1],]
rc$vote.data[order(discrimination, decreasing=TRUE)[2],]

# top 2 most "discriminatory" bills for NEGATIVE values of scale
rc$vote.data[order(discrimination)[1],]
rc$vote.data[order(discrimination)[2],]

################################################################
## COMPARING IRT WITH WNOMINATE
################################################################

library(wnominate)
nom <- wnominate(rc, dims=2, polarity=c("Cruz (R-TX)", "Cochran (R-MS)"))

par(mfrow=c(1,2))
plot(irt$xbar[,1], nom$legislators$coord1D, 
    xlab="IRT ideal point (1D)", ylab="W-NOMINATE (1D)")
plot(irt$xbar[,2], nom$legislators$coord2D, 
    xlab="IRT ideal point (2D)", ylab="W-NOMINATE (2D)")


################################################################
## IRT MODEL FIT
################################################################

## 1) PROPORTION OF CORRECTLY PREDICTED VOTES

# baseline
repub <- ifelse(rc$legis.data$party=="R", 0, 1) # dummy 'legislator==Republican'
K <- dim(rc$votes)[2] # number of votes
tab <- table(c(rc$votes), rep(repub, K)) # baseline: all Rs vote together
sum(diag(tab))/sum(tab) # proportion of correctly predicted

# W-NOMINATE
nom$fits[1:2]

# 1-dimensional model
pred <- matrix(NA, nrow=dim(irt$xbar)[1], ncol=dim(irt$betabar)[1]) # empty matrix
for (i in 1:nrow(pred)){
  for (j in 1:ncol(pred)){
    # compute predicted probability that legislator i votes YES to bill j
    pred[i,j] <- plogis(irt$xbar[i,1] * irt$betabar[j,1] - irt$betabar[j,3])
  }
}
tab <- table(c(rc$votes), c(pred)>0.50)
sum(diag(tab))/sum(tab)

# 2-dimensional model
pred <- matrix(NA, nrow=dim(irt$xbar)[1], ncol=dim(irt$betabar)[1])
for (i in 1:nrow(pred)){
  for (j in 1:ncol(pred)){
    pred[i,j] <- plogis(irt$xbar[i,1] * irt$betabar[j,1] + 
        irt$xbar[i,2] * irt$betabar[j,2] - irt$betabar[j,3])
  }
}
tab <- table(c(rc$votes), c(pred)>0.50)
sum(diag(tab))/sum(tab)

## 2) PROPORTION OF 'YES' VOTES BY PROBABILITY BINS
# (useful for sparse vote matrices)
bins <- mapply(function(x, y) which(c(pred)>x & c(pred)<=y), 
            seq(0, .9, .1), seq(.1, 1, .1))
pred.bins <- lapply(bins, function(x) mean(c(rc$votes)[x]==1, na.rm=TRUE))

plot(seq(.05, .95, .10), pred.bins, xlab="Probability bins", ylab="% Yeas")
lines(seq(.05, .95, .10), pred.bins)
abline(a=0, b=1)


## 3) PROPORTION OF CORRECTLY PREDICTED VOTES USING ESTIMATED CUTPOINTS

# function to compute correctly predicted votes for a single bill
# for a given cutpoint
max.pred <- function(cutpoint, vote, xbar){
    tab <- table(xbar > cutpoint, vote)
    pred <- sum(diag(tab)) / sum(tab)
    ifelse(pred>0.50, pred, 1-pred)
}

# example
max.pred(-1, rc$votes[,1], irt$xbar[,1])

# function to compute cutpoints and % correctly predicted votes
compute.cutpoints <- function(ideal.points, votes){
    # loop over votes
    cutpoints <- apply(votes, 2, function(x)
        optimize(max.pred, interval=range(ideal.points, na.rm=TRUE),
            vote=x, xbar=ideal.points, maximum=TRUE))
    return(matrix(unlist(cutpoints), ncol=2, byrow=TRUE))
}

cutpoints <- compute.cutpoints(irt$xbar[,1], rc$votes)
mean(cutpoints[,2])

cutpoints <- compute.cutpoints(irt$xbar[,2], rc$votes)
mean(cutpoints[,2])


################################################################
## IRT WITH STAN
################################################################

library(rstan)

stan.code <- '
data {
  int<lower=1> J; // number of legislators
  int<lower=1> K; // number of bills
  int<lower=1> N; // number of observations
  int<lower=1,upper=J> j[N]; // legislator for observation n
  int<lower=1,upper=K> k[N]; // bill for observation n
  int<lower=0,upper=1> y[N]; // vote of observation n
}
parameters {
  real alpha[K];             
  real beta[K];   
  real theta[J];
}
model {
  alpha ~ normal(0, 25);
  beta ~ normal(0, 25);
  theta ~ normal(0, 1);
  for (n in 1:N)
    y[n] ~  bernoulli_logit( theta[j[n]] * beta[k[n]] - alpha[k[n]] );
}
'

J <- dim(rc$votes)[1]
K <- dim(rc$votes)[2]
N <- length(rc$votes)
j <- rep(1:J, times=K)
k <- rep(1:K, each=J) 
y <- c(rc$votes)

# deleting missing values
miss <- which(is.na(y))
N <- N - length(miss)
j <- j[-miss]
k <- k[-miss]
y <- y[-miss]

## data and initial values
stan.data <- list(J=J, K=K, N=N, j=j, k=k, y=y)
inits <- list(list(alpha=rnorm(K, 0, 2), beta=rnorm(K, 0, 2),
    theta=ifelse(rc$legis.data$party=="R", 1, -1)))

stan.fit <- stan(model_code=stan.code, data=stan.data, iter=500, warmup=200,
    chains=1, thin=2, inits=inits)

load("lab14_stan_irt.Rdata")


## convergence
traceplot(stan.fit, pars='theta', ask=TRUE)

## comparing with WNOMINATE and Jackman's ideal
estimates <- summary(stan.fit)
theta <- estimates$summary[paste0("theta[", 1:J, "]"),1]

par(mfrow=c(1,2))
plot(irt$xbar[,1], theta,
    xlab="IRT ideal point (1D)", ylab="IRT ideal point (STAN)")
plot(nom$legislators$coord1D, theta,
    xlab="W-NOMINATE (1D)", ylab="IRT ideal point (STAN)")


################################################################
## IRT WITH COVARIATES
################################################################

# Different ways of doing this...
# With Stan, it would be something like this:

stan.code <- '
data {
  int<lower=1> J; // number of legislators
  int<lower=1> K; // number of bills
  int<lower=1> N; // number of observations
  int<lower=1,upper=J> j[N]; // legislator for observation n
  int<lower=1,upper=K> k[N]; // bill for observation n
  int<lower=0,upper=1> y[N]; // vote of observation n
  real party[J]; // party of legislator j (0 for D/I, 1 for R)
}
parameters {
  real alpha[K];             
  real beta[K];   
  real theta[J]; # realized ideology
  real gamma[J]; # unobserved, true ideology
  real beta_party; # effect of party ID
  real<lower=0.1> sigma; # sd of legislator ideology
}
model {
  alpha ~ normal(0, 5); 
  beta ~ normal(0, 5);
  beta_party ~ normal(0, 2);
  for (i in 1:J){
    theta[i] ~ normal(gamma[i], sigma); # hierarchical structure for obs. ideol.
    gamma[i] ~ normal(beta_party * party[i], sigma); # effect of party ID
  };
  for (n in 1:N)
    y[n] ~  bernoulli_logit( theta[j[n]] * beta[k[n]] - alpha[k[n]] );
}
'

repub <- ifelse(rc$legis.data$party=="R", 0, 1)
stan.data <- list(J=J, K=K, N=N, j=j, k=k, y=y, party=repub)

inits <- list(list(alpha=rnorm(K, 0, 2), beta=rnorm(K, 0, 2),
    theta=ifelse(rc$legis.data$party=="R", 1, -1),
    gamma=ifelse(rc$legis.data$party=="R", 1, -1),
    beta_party=1, sigma=0.5))

stan.fit <- stan(model_code=stan.code, data=stan.data, iter=500, warmup=200,
    chains=1, thin=2, inits=inits)


# MCMCpack has a function to do this too
female <- ifelse(rc$legis.data$gender=="F", 1, 0)
mcmc <- MCMCirtHier1d(rc$votes, data.frame(republican=repub, female=female))

results <- summary(mcmc)
round(results$statistics['beta.republican',],2)
round(results$statistics['beta.female',], 2)
theta <- results$statistics[1:104,1]

plot(irt$xbar[,1], theta,
    xlab="IRT ideal point (1D)", ylab="IRT ideal point with covariates (MCMCpack)")


