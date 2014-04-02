################################################################
## Bayesian Statistics: Fitting a Probit Regression
## Quant III Lab 13
## December 5th 2013
################################################################

# Simon Jackman's R package
#install.packages("pscl")

# Additional functions to monitor convergence of chains 
#install.packages("coda")

# Function to convert stan fit objects to coda objects
stan_to_coda <- function(fit){
    t <- extract(fit, permuted=FALSE, inc_warmup=FALSE)
    mcmc <- mcmc.list(lapply(1:ncol(t), function(x) mcmc(t[,x,])))
    return(mcmc)
}

library(pscl)
library(coda)

data(iraqVote)
str(iraqVote)

# example in page 391 of Jackman

# codebook:
?iraqVote

# Frequentist probit regression using MLE
summary(glm(y ~ gorevote + rep, data=iraqVote, family = binomial("probit")))

# Bayesian regression with Stan
library(rstan)

probit_code <- '
    data {
        int<lower=0> N; # observations
        int<lower=0> K; # variables
        row_vector[K] x[N];  # covariates
        int<lower=0,upper=1> y[N]; # outcome variable
    }
    parameters {
        vector[K] beta;  ## coefficients to be estimated
    }
    model {
        for (k in 1:K)
            beta[k] ~ normal(0,1000);  # priors
        for (n in 1:N)
            y[n] ~ bernoulli(Phi(x[n] * beta));  # distribution of outcome variable
    }
'

x <- cbind(1, iraqVote$gorevote, iraqVote$rep)
y <- iraqVote$y
data <- list(N=length(y), K=dim(x)[2], y=y, x=x)
fit <- stan(model_code=probit_code, data=data, iter=1000, chains=4)

# summary results
monitor(fit, digits_summary=2)

# Assessing convergence
# 1) Looking at traceplots
traceplot(fit, pars='beta')

# 2/3) effective sample size and R hats
monitor(fit, digits_summary=2)

# 4) Geweke tests
mcmc <- stan_to_coda(fit)
geweke.diag(mcmc)
geweke.plot(mcmc, ask=TRUE)

# 5) Heidelberger-Welch test of non-stationarity
heidel.diag(mcmc)

# 6) looking at serial correlation in chains
autocorr(mcmc)

autocorr.plot(mcmc, ask=TRUE)
