################################################################
## Time Series and Stationarity
## Author: Pablo Barberá
## Quant III Lab 8
## October 24th 2013
################################################################

#install.packages("urca")

## Let's start with the following model in mind:
## y_t = rho * y_{t-1} + nu
## where nu is iid

## now, if |rho| < 1, then we have a stationary time-series
## if |rho|=1, then we have a non-stationary time-series (random walk)
## if |rho|>1, then we have a non-stationary time-series (explosive)

## today we will focus on the second situation, also known as unit root

## remember from the lectures that in this type of time series the current
## value of y is simply the sum of all past random shocks. As a result,
## the series tends to drift up and down (although it ISN'T trending), 
## although the series is mean-reverting in the long run (it's just the sum
## of iid errors)

## However, the variance of a random walk does increase with time

## This R script illustrates these points and shows how to test for
## unit roots

################################################################
## y is I(1), x is stationary, unit root
################################################################

set.seed(123)
t <- 500
x <- 5 + rnorm(t)
nu <- rnorm(t)
e <- rnorm(t)
rho <- 1 ## unit root
for (i in 2:t){ e[i] <- rho * e[i-1] + nu[i] }
beta <- 0.25
y <- beta * x + e

par(mfrow=c(1,2))
## y is non-stationary time series (random walk)
plot(y, type="l")
## x is stationary
plot(x, type="l")

## inconsistent estimates
summary(lm(y ~ x))

## dickey-fuller test for unit root
library(urca)
summary(ur.df(y, type="none", lags=0))
# value of test-statistic is HIGHER than critical values,
# (careful with negative signs)
# so we cannot reject null hypothesis of unit root
# i.e., this test SUGGESTS we might have unit root
# (see chapter 5 from "ts_part.pdf", page 114)

## in STATA, you would type "dfuller y, noconstant"

## check also the plot that ur.df() generates:
plot(ur.df(y, type="none", lags=0))

## there are many flavors of the DF test:
## with drift / with trend / with higher order lags ("augmented DF test")

## integrating the series
d.y <- diff(y)
plot(1:(t-1), d.y, type="l")

## now we get consistent testimates
summary(lm(diff(y) ~ diff(x)))
## and no unit root
summary(ur.df(d.y, type="none", lags=0)) ## reject unit root


################################################################
## x is I(1), y is I(1), error is white noise
################################################################

set.seed(123)
t <- 500
x <- rnorm(t)
for (i in 2:t){ x[i] <- x[i-1] + rnorm(1) }
beta <- 3
y <- beta * x + rnorm(t)

## x is non-stationary (unit root)
par(mfrow=c(1,2))
summary(ur.df(x, type="none", lags=0))
plot(x, type="l")
## y is ALSO non-stationary
plot(y, type="l")
summary(ur.df(y, type="none", lags=0))
## why? the sum of an I(0) and a I(1) series is I(1)

## standard error is wrong
summary(lm(y ~ x))

## we need to co-integrate
summary(lm(diff(y) ~ diff(x)))

################################################################
## spurious regression
################################################################

## this is the example from class: two unrelated nonstationary series, if we 
## use OLS we find a strong correlation

set.seed(333)
t <- 500
x <- rnorm(t)
for (i in 2:t){ x[i] <- x[i-1] + rnorm(1) } ## x is I(1)
y <- rnorm(t)
for (i in 2:t){ y[i] <- y[i-1] + rnorm(1) } ## y is also I(1)

## both non-stationary
summary(ur.df(x, type="none", lags=0))
summary(ur.df(y, type="none", lags=0))
par(mfrow=c(1,2))
plot(y, type="l")
plot(x, type="l")

## spurious regression
summary(lm(y ~ x))







