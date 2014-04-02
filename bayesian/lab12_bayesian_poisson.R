################################################################
## Bayesian Statistics: Fitting a Bayesian Poisson Regression
## Quant III Lab 12
## November 21st 2013
################################################################

# Consider the following model: y_i ~ Poisson(exp(X_i \beta))
# for i = 1, ..., n
# Assume independent Cauchy prior distributions with location 0
# and scale 2.5 on the elements of \beta

# why Cauchy?
# http://www.stat.columbia.edu/~gelman/research/published/priors11.pdf
# equivalent to t-student distribution with df=1

x <- seq(-5, 5, .01)
par(mar=c(3,3,2,3))
plot(x, (dcauchy(x, location=0, scale=2.5)), type="l")
mtext("density", 2, line=2)

#############################################################
### Metropolis algorithm for a poisson regression ###
#############################################################

# function to compute log posterior density
log.post.dens <- function(beta, y, x){
  prior.dens <- sum(dcauchy(beta, location=0, scale=2.5, log=TRUE))
  lk <- sum(dpois(y, exp(x%*%beta), log=TRUE))
  return(prior.dens + lk) # Bayesian mantra
}

# metropolis-hastings sampler
metropolis.poisson <- function(y, x, iters=100, delta=0.25, chains=4){

  # preparing vector for stored samples
  beta.samples <- array(NA, dim=c(iters, chains, dim(x)[2]))
  lpd <- array(NA, dim=c(iters, chains))

  # chains of the metropolis algorithm
  for (chain in 1:chains){
    # drawing starting points
    beta.cur <- rnorm(n=dim(x)[2], mean=0, sd=1)
    
    # loop over iterations
    for (iter in 1:iters){
      # sampling proposal values
      beta.cand <- sapply(beta.cur, function(x) runif(n=1, min=x-delta, max=x+delta))
      # computing acceptance probability
      accept.prob <- exp(log.post.dens(beta.cand, y, x) - log.post.dens(beta.cur, y, x))
      alpha <- min(accept.prob, 1)
      # jumping with probability alpha
      if (runif(1)<=alpha) { beta.cur <- beta.cand}
      # storing samples
      beta.samples[iter, chain,] <- beta.cur
      lpd[iter,chain] <- log.post.dens(beta.cur, y, x)
    }

  }
  return(list(samples=beta.samples, lpd=lpd))
}

#############################################################
### Running the Metropolis Algorithm                      ###
#############################################################

set.seed(777)

# simulating data
library(MASS)
x <- mvrnorm(n=50, mu=c(-0.5,0,0.5), Sigma=diag(3))
true.beta <- c(-1, 0, 1)
y <- rpois(n=50, lambda=exp(x%*%true.beta))

# fitting model
fit <- metropolis.poisson(y, x, iters=500, delta=0.10, chains=4)

# checking convergence of chains and getting summary of results
library(rstan)
monitor(fit$samples)

# trace plots
par(mfrow=c(2,2), mar=c(3, 3, 2, 3))
colors <- c("black", "red", "blue", "green")
for (par in 1:3){
  plot(1:500, fit$samples[1:500,1,par], type="l", ylim=c(-1.5,1.5))
  for (chain in 2:4){
    lines(1:500, fit$samples[1:500,chain,par], type="l", col=colors[chain])
  }
mtext("Iteration", side=1, line=2)
mtext(paste("beta[", par, "]", sep=""), side=2, line=2)
}
plot(1:500, fit$lpd[1:500,1], type="l")
  for (chain in 2:4){
    lines(1:500, fit$lpd[1:500,chain], type="l", col=colors[chain])
  }
mtext("Iteration", side=1, line=2)
mtext("Log Posterior Density", side=2, line=2)


# checking with frequentist poisson regression

summary(glm.fit <- glm(y ~ -1 + x, family=poisson(link=log) ))



#############################################################
### Hamiltonian Monte-Carlo                               ###
#############################################################


library(rstan)

hmc_code <- '
    data {
        int<lower=0> N; # observations
        int<lower=0> K; # variables
        int y[N];       # data (integers)
        row_vector[K] x[N];  # covariates
    }
    parameters {
        vector[K] beta;  ## coefficients to be estimated
    }
    model {
      for (k in 1:K)
        beta[k] ~ cauchy(0, 2.5); ## priors on betas
      for (n in 1:N)
        y[n] ~ poisson(exp(x[n] * beta)); ## model for y
    }
'

data <- list(N=length(y), K=dim(x)[2], y=y, x=x)
fit <- stan(model_code=hmc_code, data=data, iter=500, chains=4)

# summary results
monitor(fit)

# convergence
traceplot(fit, pars='beta')

# what happens when you change priors to, for example, normal(0, 0.05)



#############################################################
### Another example with HMC                              ###
#############################################################

# Bayesian regression (from Stan manual, page 87)

x <- runif(100)
beta <- 3
alpha <- 1
y <- alpha + beta * x + rnorm(100)

hmc_code <- '
    data {
        int<lower=0> N; # observations
        vector[N] x; # covariate
        vector[N] y; # outcome variable
    }
    parameters {
        real alpha;  # intercept
        real beta;   # slope
        real<lower=0> sigma; # error
    }
    model {
        alpha ~ normal(0, 100);
        beta ~ normal(0, 100);
        for (n in 1:N)
          y[n] ~ normal(alpha + beta * x[n], sigma);
    }
'

data <- list(N=length(y), x=x, y=y)
fit <- stan(model_code=hmc_code, data=data, iter=500, chains=4)

monitor(fit)



