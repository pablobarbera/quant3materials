################################################################
## Bayesian Statistics: Fitting a Hierarchical Regression Model
## Quant III Lab 13
## December 5th 2013
################################################################

# Predicting baseball batting averages via a hierarchical model
# (classic examples that shows Bayesian models outperforming MLE)
# See pages 310-316 Jackman

library(pscl)
library(rstan)

# Function to convert stan fit objects to coda objects
stan_to_coda <- function(fit){
    t <- extract(fit, permuted=FALSE, inc_warmup=FALSE)
    mcmc <- mcmc.list(lapply(1:ncol(t), function(x) mcmc(t[,x,])))
    return(mcmc)
}

data(EfronMorris)

d <- EfronMorris

# codebook:
?EfronMorris

# MLE estimate is just batting average
rmse <- sqrt(mean((d$y-d$p)^2))
rmse

# Hierarchical model ('partial pooling')

# each observed batting average y_i is drawn from a normal density with
# unknown mean theta_i; and fixed variance sigma:
sigma <- sqrt(mean(d$y) * (1-(mean(d$y)))/45)

# Each theta_i (the true batting average of each player; which we would observe
# if they kept playing fover) is also drawn from a normal distribution with
# unknown mean mu and variance tau.

# Then we give priors to the hyperparameters mu and tau.
# for mu, we know it must be somewhere between 0.15 and .30; so we choose values
# that give a 95% credible interval: normal(0.225, 0.0375)
x <- seq(0, 0.5, 0.001)
plot(x, dnorm(x, mean=0.225, sd=0.0375), type="l")

# for tau, we choose a gamma distribution (why? it's always positive
# and it's a conjugate prior for a normal distribution)

# what prior belief do we have about tau, the standard deviation of the thetas?
# the sd of the observed y is:
sd(d$y)

# it must be somewhere around 0.07; probably higher; let's say 0.10
# so MEAN of values we draw from gamma distribution must be around 0.10
# gamma distribution has mean shape/rate [OR shape*scale]

# so let's follow Jackman and choose shape=7; then rate=7/0.10 ~= 70
x <- seq(0.001, 1, 0.001)
plot(x, dgamma(x, shape=7, rate=70), type="l")

hierarchical_code <- '
    data {
        int<lower=0> N; # observations
        real y[N]; # outcome variable
        real<lower=0> sigma; # fixed variance of observed batting averages
    }
    parameters {
        real theta[N]; ## unobserved true batting averages
        real mu; ## hyperparameter: mean of batting averages
        real<lower=0> tau; ## hyperparameter: variance of batting averages
    }
    model {
        mu ~ normal(0.225, 0.0375); ## prior about hyperparameter (overall batting average)
        tau ~ gamma(7, 70); ## priors about variance of overall batting average
        for (n in 1:N){;
            y[n] ~ normal(theta[n], sigma);
            theta[n] ~ normal(mu, tau);
        }
    }
'

data <- list(N=length(d$y), y=d$y, sigma=sigma)
fit <- stan(model_code=hierarchical_code, data=data, iter=1000, chains=4)

# summary results
monitor(fit, digits_summary=3)

# Assessing convergence
# 1) Looking at traceplots
traceplot(fit, pars='theta', inc_warmup=FALSE, ask=TRUE)

# 4) Geweke tests
mcmc <- stan_to_coda(fit)
geweke.diag(mcmc)

# 5) Heidelberger-Welch test of non-stationarity
heidel.diag(mcmc)

# 6) looking at serial correlation in chains
autocorr(mcmc)

autocorr.plot(mcmc, ask=TRUE)


# computing rmse
estimates <- summary(fit)
predicted <- estimates$summary[1:18, 'mean']
rmse <- sqrt(mean((predicted-d$p)^2))
rmse

# visualizing results (replication of Figure 7.4 in page 316 Jackman)
actual <- data.frame(batter = d$name, average=d$p, estimate="Actual")
bayes <- data.frame(batter = d$name, average=predicted, estimate="Bayes")
mle <- data.frame(batter = d$name, average=d$y, estimate="MLE")

df <- rbind(actual, bayes, mle)

library(ggplot2)
p <- ggplot(df, aes(x=estimate, y=average, group=batter))
p + geom_point() + geom_line() + theme_bw() + coord_flip() +
    scale_x_discrete(expand=c(0.05,0.05)) +
    scale_y_continuous("Batting average") + 
    theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(size = 0.50)) +
    geom_vline(xintercept=1, alpha=1/3, size=0.2) +
    geom_vline(xintercept=2, alpha=1/3, size=0.2) +
    geom_vline(xintercept=3, alpha=1/3, size=0.2)







