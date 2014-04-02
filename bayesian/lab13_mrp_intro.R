################################################################
## Bayesian Statistics: Fitting a Hierarchical Regression Model
## (Basic example of multilevel regression with poststratification)
## Quant III Lab 13
## December 5th 2013
################################################################

# Inspired by: http://www.princeton.edu/~jkastell/mrp_primer.html

library(foreign)
library(rstan)

## polls about gay rights ("ground truth"), from Lax and Phillips (2009)
d <- read.dta("lax_phillips_gay_rights.dta", convert.underscore = TRUE) 

# delete missing data in variables of interest
d <- d[!is.na(d$yes.of.all) & !is.na(d$age) & (d$state!="") & !is.na(d$region),]

# computing means by state
actual.means <- aggregate(d$yes.of.all, by=list(d$state), mean)$x

# computing average age by state
age_means <- aggregate(d$age, by=list(d$state), mean)$x

# take random sample of N=1000 (for computation purposes)
set.seed(1234)
d <- d[sample(1:length(d$state), 1000),]

# and let's see how well we can replicate "true" support for gay marriage
# in each state using a random sample

# Variables at the individual level
y <- d$yes.of.all 
age <- d$age
state <- d$state
states <- sort(unique(state))
ss <- match(state, states) ## number of state to which each obs. belongs

################################################################
## RANDOM INTERCEPTS MODEL
################################################################

multilevel_code <- '
    data {
        int<lower=0> N; # number of respondents
        int<lower=0> S; # number of states
        int<lower=0, upper=1> y[N]; # outcome variable
        int<lower=1,upper=S> ss[N]; # state for each observation
        real age[N];
        real age_means[S];
    }
    parameters {
        real alpha[S]; # state-specific intercept
        real beta; # effect of age
        real mu; # mean of intercepts
        real<lower=0> tau; # sd of intercepts
    }
    model {
        alpha ~ normal(0, 100);
        beta ~ normal(0, 100);
        mu ~ normal(0, 100);
        tau ~ normal(0, 100);
        for (s in 1:S)
            alpha[s] ~ normal(mu, tau);
        for (n in 1:N)
            y[n] ~ bernoulli(inv_logit(alpha[ss[n]] + beta * age[n]));
    }
    generated quantities {
        real predictions[S];
        for (s in 1:S)
            predictions[s] <- inv_logit(alpha[s] + beta * age_means[s]);
    }
'
data <- list(N=length(y), S=length(states), y=y,
    ss=ss, age=age, age_means=age_means)

fit <- stan(model_code=multilevel_code, data=data, iter=1000, chains=2)

# backup
load("lab13_irt_results.rda")

# results and convergence statistics
monitor(fit, digits_summary=2, warmup=100)

traceplot(fit, ask=TRUE)

# compute observed means in random sample of N=1000
observed.means <- aggregate(y, by=list(state), mean)$x

# compute predicted means by Bayesian model
estimates <- summary(fit)
varnames <- paste0('predictions[', 1:49, ']')
predicted.means <- estimates$summary[varnames, 'mean']

# visualizing results
par(mfrow=c(1,2), mar=c(3, 3, 2, 3))
plot(observed.means, actual.means, type="n")
text(observed.means, actual.means, labels=states)

plot(predicted.means, actual.means, type="n")
text(predicted.means, actual.means, labels=states)

cor(observed.means, actual.means)
sqrt(mean((observed.means-actual.means)^2))

cor(predicted.means, actual.means)
sqrt(mean((predicted.means-actual.means)^2))

################################################################
## MODELING RANDOM INTERCEPTS
################################################################

# Variable at the state level: vote for Kerry in 2004 in each state
kerry <- c(36.7999992370605, 44.5999984741211, 44.4000015258789, 54.2999992370605, 
47, 54.2999992370605, 89.1999969482422, 53.4000015258789, 47.0999984741211, 
41.4000015258789, 49.2000007629395, 30.2999992370605, 54.7999992370605, 
39.2999992370605, 36.5999984741211, 39.7000007629395, 42.2000007629395, 
61.9000015258789, 55.9000015258789, 53.5999984741211, 51.2000007629395, 
51.0999984741211, 46.0999984741211, 39.7999992370605, 38.5999984741211, 
43.5999984741211, 35.5, 32.7000007629395, 50.2000007629395, 52.9000015258789, 
49.0999984741211, 47.9000015258789, 58.4000015258789, 48.7000007629395, 
34.4000015258789, 51.4000015258789, 50.9000015258789, 59.4000015258789, 
40.9000015258789, 38.4000015258789, 42.5, 38.2000007629395, 26, 
45.5, 58.9000015258789, 52.7999992370605, 49.7000007629395, 43.2000007629395, 
29.1000003814697)

multilevel_code <- '
    data {
        int<lower=0> N; # number of respondents
        int<lower=0> S; # number of states
        int<lower=0, upper=1> y[N]; # outcome variable
        int<lower=1,upper=S> ss[N]; # state for each observation
        real age[N];
        real kerry[S];
        real age_means[S];
    }
    parameters {
        real alpha[S]; # state-specific intercept
        real beta; # effect of age
        real gamma; # effect of kerry
        real mu; # mean of intercepts
        real<lower=0> tau; # sd of intercepts
    }
    model {
        alpha ~ normal(0, 100);
        beta ~ normal(0, 100);
        gamma ~ normal(0, 100);
        mu ~ normal(0, 100);
        tau ~ normal(0, 100);
        for (s in 1:S)
            alpha[s] ~ normal(mu + gamma * kerry[s], tau);
        for (n in 1:N)
            y[n] ~ bernoulli(inv_logit(alpha[ss[n]] + beta * age[n]));
    }
    generated quantities {
        real predictions[S];
        for (s in 1:S)
            predictions[s] <- inv_logit(alpha[s] + beta * age_means[s]);
    }
'


data <- list(N=length(y), S=length(states), y=y,
    ss=ss, age=age, age_means=age_means, kerry=kerry)

fit <- stan(model_code=multilevel_code, data=data, iter=1000, chains=2)

# results and convergence statistics
monitor(fit, digits_summary=2, warmup=100)
#traceplot(fit, ask=TRUE)

# compute predicted means by Bayesian model
estimates <- summary(fit)
varnames <- paste0('predictions[', 1:49, ']')
predicted.means <- estimates$summary[varnames, 'mean']

# visualizing results
par(mfrow=c(1,2), mar=c(3, 3, 2, 3))
plot(observed.means, actual.means, type="n")
text(observed.means, actual.means, labels=states)

plot(predicted.means, actual.means, type="n")
text(predicted.means, actual.means, labels=states)

cor(observed.means, actual.means)
sqrt(mean((observed.means-actual.means)^2))

cor(predicted.means, actual.means)
sqrt(mean((predicted.means-actual.means)^2))




