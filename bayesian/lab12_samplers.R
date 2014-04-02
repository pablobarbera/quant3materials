################################################################
## Bayesian Statistics: Sampling algorithms
## Quant III Lab 12
## November 21st 2013
################################################################

# install.packages("mvtnorm")
library(mvtnorm)
library(rstan)

# Suppose a single observation (y1, y2) from a bivariate normally distributed
# population with unknown mean theta = (theta1, theta2) and known covariance
# matrix (1, rho / rho, 1) -- so covariance between theta1 and theta2 is rho

# Assume a uniform prior distribution on theta.

# How can we sample from the posterior distribution of theta | y ?

# Posterior distribution is:
# (theta1, theta2) | y ~ N( (y1, y2) , (1, rho / rho, 1))

# This example comes from Gelman, BDA 3, page 277
# You can find a similar example in Jackman, page 214-215
# Also in Chib and Greenberg, 1995, page 333

# Data
y <- c(0, 0)
rho <- 0.9
Sigma <- matrix(c(1, rho, rho, 1), nrow=2)

# And note that this is a trivial example! We could just use one of R's
# built-in functions (it's also easy to do manually in the Cholesky approach)

theta <- rmvnorm(1000, mean=y, sigma=Sigma)
plot(theta[,1], theta[,2], cex=.50, pch=20)
mtext("theta1", side=1, line=1.75)
mtext("theta2", side=2, line=1.75)

# But in general they are not available, usually we just have the density
# function of our target distribution and the data, and we have to work with
# that, either programming a sampler or using Stan.


################################################################
### GRID SAMPLING
################################################################

# This is a 'brute force' approach: we just take a grid of possible values,
# compute the posterior density and sample with probability equal to that
# posterior density.

# simulated parameter grid
n.grid <- 100
theta1 <- seq(-3, 3, length.out=n.grid)
theta2 <- seq(-3, 3, length.out=n.grid)

# width of space between possible values of the parameters
width <- theta1[2] - theta1[1]

library(mvtnorm)

# 'fill' grid of possible parameter values with posterior density
posterior <- matrix(NA, nrow=n.grid, ncol=n.grid)

for (i in 1:n.grid){
    for (j in 1:n.grid){
        posterior[i, j] <- 1 * dmvnorm(y, mean=c(theta1[i], theta2[j]), sigma=Sigma)
    }
}

# (note that I multiply by 1 simply to illustrate that if you would use priors,
# there you would have the prior probability. In this case, we use uniform
# priors, so any constant will give us the same answer.)

# normalizing (so that max of posterior is 1)
posterior <- posterior / max(posterior)

# contour plot
contours <- c(seq(.05, .95, .10), .99)
par(mar=c(3, 3, 2, 3))
contour (theta1, theta2, posterior, levels=contours,
cex=2, drawlabels=TRUE)
mtext("theta1", side=1, line=1.75)
mtext("theta2", side=2, line=1.75)

# computing marginal distribution of alpha: p(theta1|theta2,y)
post.theta1 <- rowSums(posterior) / sum(rowSums(posterior))

# sample from conditional distribution
n.sim <- 1000
theta1.sampled <- sample(theta1, size=n.sim, prob=post.theta1, replace=TRUE)

# draw theta2 from p(theta2|theta1,y)
theta2.sampled <- rep(NA, n.sim)
for (i in 1:n.sim){
    loc <- which(theta1==theta1.sampled[i])
    theta2.sampled[i] <- sample(theta2, size=1, prob=posterior[loc,])
}

# plot samples (with jitter -- otherwise we would only see values in grid!)
theta1.sampled <- theta1.sampled + runif(n.sim, min=0, max=width)
theta2.sampled <- theta2.sampled + runif(n.sim, min=0, max=width)

# plotting draws from posterior distribution
plot(theta1.sampled, theta2.sampled, cex=.50, pch=20)
mtext("theta1", side=1, line=1.75)
mtext("theta2", side=2, line=1.75)


#### GRID SAMPLING: SUMMARY ####
# PROS: very easy for simple examples; no convergence issues.
# CONS: extremely inefficient for anything beyond 2 parameters; and also
# very imprecise unless width of grid spaces is very small.


################################################################
### METROPOLIS ALGORITHM (with random walk)
################################################################

# This was the original MCMC algorithm. Parameters start with some (usually
# random) values. From there, jump somewhere else (with a given jumping
# distribution; if uniform from 0 to delta then it's called 'random walk'). 
# Compute posterior density for new parameter values ('candidate' values). 
# Accept new position with probability r = lpd.candidate / lpd.previous. 
# If accept, then move to those values and repeat. If not, stay where you 
# are and try jumping again.


## FUNCTIONS
log.post.dens <- function(theta, y, Sigma){
    prior <- 1 # uniform prior
    likelihood <- dmvnorm(y, mean=theta, sigma=Sigma)
    return(log(prior * likelihood)) # Bayesian mantra
}

metropolis.bivariate <- function(y, Sigma, iters=100, delta=0.10, chains=4){

    # preparing vector for stored samples
    theta.samples <- array(NA, dim=c(iters, chains, length(y)))
    dimnames(theta.samples)[[3]] <- paste0('theta[', 1:length(y), ']')

    # chains of the metropolis algorithm
    for (chain in 1:chains){
    # drawing starting points
    theta.cur <- rnorm(n=length(y), mean=0, sd=1)
    
    # loop over iterations
    for (iter in 1:iters){
      # sampling proposal values
      theta.cand <- sapply(theta.cur, function(x) runif(n=1, min=x-delta, max=x+delta))
      # computing acceptance probability
      accept.prob <- exp(log.post.dens(theta.cand, y, Sigma) - log.post.dens(theta.cur, y, Sigma))
      alpha <- min(accept.prob, 1)
      # jumping with probability alpha
      if (runif(1)<=alpha) { theta.cur <- theta.cand}
      # storing samples
      theta.samples[iter, chain,] <- theta.cur
    }
  }
  return(theta.samples)
}

## RUNNING THE SAMPLER
iters <- 500
metr <- metropolis.bivariate(y, Sigma, iters=iters, delta=0.10, chains=4)

# summary of results
monitor(metr) ## it didn't converge!

# trace plots
par(mfrow=c(1,2), mar=c(3, 3, 2, 3))
colors <- c("black", "red", "blue", "green")
for (par in 1:2){
  plot(1:iters, metr[1:iters,1,par], type="l", ylim=range(metr))
  for (chain in 2:4){
    lines(1:iters, metr[1:iters,chain,par], type="l", col=colors[chain])
  }
mtext("Iteration", side=1, line=2)
mtext(paste("theta[", par, "]", sep=""), side=2, line=2)
}

# visualizing sequence of Metropolis sampler
plot(metr[1,1,'theta[1]'], metr[1,1,'theta[2]'], 
    xlim=range(metr[,,'theta[1]']),
    ylim=range(metr[,,'theta[2]']),
    xlab="theta[1]", ylab="theta[2]",
    pch = 19, cex=0.50)
for (i in 1:iters){
    lines(metr[i:(i+1),1,'theta[1]'], metr[i:(i+1),1,'theta[2]'])
    message <- paste0("Iteration ", i, ". Press enter to continue.")
    invisible(readline(message))
    i <- i + 1
}

# visualizing samples (first chain)
plot(metr[,1,'theta[1]'], metr[,1,'theta[2]'], pch=19, cex=0.30)

# How can we solve it?
# 1) Play around with 'delta' parameter (maximum jump)
# 2) Increase number of iterations

#### METROPLIS SAMPLER: SUMMARY ####
# PROS: very easy to program. It works even for relatively complex densities
# CONS: it can be VERY inefficient, requires some tuning (delta parameter).



################################################################
### GIBBS SAMPLER
################################################################

# Key idea: parameter vector theta is divided into different components
# (theta1, theta2), and each iteration cycles through each of these components,
# drawing each subset conditional on the others:
# Step 1: sample (theta_1 | theta_2, ..., theta_d, y)
# Step 2: sample (theta_2 | theta_1, ..., theta_d, y)
# ...
# Step d: sample (theta_d | theta_1, theta_2, ..., theta_d-1, y)
# And repeat for T iterations.

# This works because usually it's easier to sample from conditional posterior
# distributions of parameters.

# In our example, sampling from the bivariate normal becomes a problem of
# sampling from conditional densities, that are univariate normal (easier to
# sample)

# From BDA3, appendix, page 580:
# 'The conditional distribution of any subvector of theta given the remaining
# elements is again multivariate normal. If we partition theta into subvectors
# theta = (U, V), then p(U|V) is normal with:
# E(U|V) = E(U) + cov(U,V)var(V)^{-1}(V-E(V))
# var(U|V) = var(U) - cov(U,V)var(V)^{-1}cov(V,U)

# In our example, the conditional posterior densities become:

# theta1 | theta2, y ~ N(mu, Sigma)
# with:
# mu = y[1] + rho * (theta2 - y[2])
# Sigma = 1 - rho^2

# That's what I do below

## FUNCTIONS
update.theta.1 <- function(theta1, theta2, y, rho){
    rnorm(n=1, mean=y[1]+rho*(theta2-y[2]), sd=sqrt(1-rho^2))
}

update.theta.2 <- function(theta1, theta2, y, rho){
    rnorm(n=1, mean=y[2]+rho*(theta1-y[1]), sd=sqrt(1-rho^2))
}

gibbs.bivariate <- function(y, rho, iters=100, chains=4){

    # preparing vector for stored samples
    theta.samples <- array(NA, dim=c(iters, chains, length(y)))
    dimnames(theta.samples)[[3]] <- paste0('theta[', 1:length(y), ']')

    # chains of the gibbs sampler
    for (chain in 1:chains){
        # starting points
        theta <- rnorm(n=length(y))
        for (iter in 1:iters){
            # sampling theta 1 from conditional distribution
            theta[1] <- update.theta.1(theta[1], theta[2], y, rho)
            # sampling theta 2 from conditional distribution
            theta[2] <- update.theta.2(theta[1], theta[2], y, rho)
            # storing samples
            theta.samples[iter, chain, ] <- theta
        }
    }
    return(theta.samples)
}


## RUNNING THE SAMPLER
iters <- 500

gibbs <- gibbs.bivariate(y, rho, iters=iters, chains=4)

# summary of results
monitor(gibbs) # it converged!

# trace plots
par(mfrow=c(1,2), mar=c(3, 3, 2, 3))
colors <- c("black", "red", "blue", "green")
for (par in 1:2){
  plot(1:iters, gibbs[1:iters,1,par], type="l", ylim=range(gibbs))
  for (chain in 2:4){
    lines(1:iters, gibbs[1:iters,chain,par], type="l", col=colors[chain])
  }
mtext("Iteration", side=1, line=2)
mtext(paste("theta[", par, "]", sep=""), side=2, line=2)
}

# visualizing sequence of Gibbs sampler
plot(gibbs[1,1,'theta[1]'], gibbs[1,1,'theta[2]'], 
    xlim=range(gibbs[,,'theta[1]']),
    ylim=range(gibbs[,,'theta[2]']),
    xlab="theta[1]", ylab="theta[2]",
    pch = 19, cex=0.50)
for (i in 1:iters){
    message <- paste0("Iteration ", i, ". Press enter to continue.")
    lines(gibbs[i:(i+1),1,'theta[1]'], rep(gibbs[i,1,'theta[2]'], 2))
    lines(rep(gibbs[(i+1),1,'theta[1]'], 2), gibbs[i:(i+1),1,'theta[2]'])
    invisible(readline(message))
    i <- i + 1
}

# visualizing samples (first chain)
plot(gibbs[,1,'theta[1]'], gibbs[,1,'theta[2]'], pch=19, cex=0.30)


#### GIBBS SAMPLER: SUMMARY ####
# PROS: easy-(ish?) to program. For some problems, VERY efficient. Nice way
# to split multidimensional problems into simpler densities.
# CONS: you need to do math (ugh) to get the conditional distributions. Not
# all densities can be split into conditionals with a nice equation.



################################################################
### HAMILTONIAN MONTE CARLO (with STAN)
################################################################

# You can read Gelman's book to understand how HMC works. It's not trivial.
# The basic idea is that jumping rules are much more efficient because they
# learn from the gradient of the log posterior density, so they know better
# where to jump to. As a result, it can be MUCH more efficient.

# Luckily, we do not need to program HMC because stan does (almost) all the work
# for us.


library(rstan)

hmc_code <- '
    data {
        vector[2] y; // data (means)
        matrix[2,2] Sigma; // covariance matrix (known)
    }
    parameters {
        vector[2] theta; // parameters of bivariate normal
    }
    model {
        y ~ multi_normal(theta, Sigma);
    }
'

data <- list(y=y, Sigma=Sigma)
fit <- stan(model_code=hmc_code, data=data, iter=500, chains=4)


# summary of results
monitor(fit)

# trace plots
traceplot(fit, pars='theta')

# visualizing sequence of HMC sampler
hmc <- extract(fit, permuted=FALSE, inc_warmup=TRUE)

plot(hmc[1,1,'theta[1]'], hmc[1,1,'theta[2]'], 
    xlim=range(hmc[,,'theta[1]']),
    ylim=range(hmc[,,'theta[2]']),
    xlab="theta[1]", ylab="theta[2]",
    pch = 19, cex=0.50)
for (i in 1:iters){
    message <- paste0("Iteration ", i, ". Press enter to continue.")
    lines(hmc[i:(i+1),1,'theta[1]'], hmc[i:(i+1),1,'theta[2]'])
    invisible(readline(message))
    i <- i + 1
}


# visualizing samples (first chain)
plot(hmc[,1,'theta[1]'], hmc[,1,'theta[2]'], pch=19, cex=0.30)

#### HAMILTONIAN MONTE CARLO: SUMMARY ####
# PROS: very easy to program (just write down the statistical model!). Very
# efficient in general (or at least as efficient as Gibbs). Works for all 
# types of problems.
# CONS: you need to learn how to use stan. Less control over sampler (but
# maybe it's for the best)











