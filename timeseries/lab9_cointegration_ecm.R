################################################################
## Cointegration and error-correction models
## Author: Pablo Barberá
## Quant III Lab 9
## October 31th 2013
################################################################

#install.packages("urca")

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

## when both x and y are I(1), then the hope is that they are cointegrated
## i.e., they are in a long run equilibrium relationship

## how do we know if that's the case? if residuals are stationary
reg <- lm(y~x)
summary(ur.df(reg$resid, type="none", lags=0))

plot(ur.df(reg$resid, type="none", lags=0))

plot(1:t, x, type="l", ylim=range(x,y))
lines(y, col="red")

## if that's the case, then \hat\beta is the "co-integrating" value

## but note that this only gives us the long run relationship! so, in the
## long run, y is always x * beta. They can "wander" a lot, but always
## together (example of drunk and his/her dog)

## but HOW FAST does the adjustment take place? when either x or y move
## outside of the equilibrium, how long does it take until they are again
## in equilibrium?

## ERROR CORRECTION MODEL(S)

## Engle-Granger two-step
## a) estimate cointegrating relationship
reg <- lm(y~x)
## b) compute errors
resid <- reg$resid
## c) now estimate short-run relationship
lm(diff(y) ~ -1 + diff(x) + resid[-t])
## the coefficient for the lagged residuals indicate how fast the
## adjustment takes place: so in one period, 99.2% of the adjustment takes
## place (that's what we would expect, given how we built the data)

## One-step method (DHSY)

## note that:
## diff.y = diff.x + l.resid
##        = diff.x + l.(y - x)
##        = diff.x + l.y - l.x
## (and ignoring coefficients for now, but it's all simple algebra)

## so we can just estimate:
lm(diff(y) ~ -1 + diff(x) + y[-t] + x[-t])
## and the estimate for l.y is same speed of adjustment
## (asymptotically, these two models are identical, after some algebra)


################################################################
## x is I(1), y is I(1), error is white noise, lower adjustment speed
################################################################

# let's simulate some data with lower adjustment speed
# i.e., it takes more than one period to "re-equilibrate"

set.seed(123)
t <- 500
x <- rnorm(t)
for (i in 2:t){ x[i] <- x[i-1] + rnorm(1) }
beta <- 3

e <- rnorm(t)
rho <- 0.7
for (i in 2:t){ e[i] <- rho * e[i-1] + e[i] }

y <- beta * x + e + 2

## y is still non-stationary (unit root)
summary(ur.df(y, type="none", lags=0))
plot(y, type="l")

## x is also non-stationary (unit root)
summary(ur.df(x, type="none", lags=0))
lines(x, col="red")

# and they are still cointegrated
reg <- lm(y~x)
summary(ur.df(reg$resid, type="none", lags=0))

# Estimating ECM:

# Engle-Granger
reg <- lm(y~x)
resid <- reg$resid
lm(diff(y) ~ -1 + diff(x) + resid[-t])
# the relationship between y and x returns to its equilibrium levels at a
# rate of about 30% of the disequilibrium each period

# DHSY
lm(diff(y) ~ -1 + diff(x) + y[-t] + x[-t])


################################################################
## x is I(1), y is I(1), error is I(1)
################################################################

# what happens when x and y are not co-integrated?

set.seed(777)
t <- 500
x <- rnorm(t)
for (i in 2:t){ x[i] <- x[i-1] + rnorm(1) } ## x is I(1)
y <- rnorm(t)
for (i in 2:t){ y[i] <- y[i-1] + rnorm(1) } ## y is also I(1)


## y is still non-stationary (unit root)
summary(ur.df(y, type="none", lags=0))
plot(y, type="l", ylim=range(x, y))

## x is also non-stationary (unit root)
summary(ur.df(x, type="none", lags=0))
lines(x, col="red")

# as we saw last week, spurious regression...
summary(lm(y~x))

# and it's not cointegrated because...
reg <- lm(y~x)
summary(ur.df(reg$resid, type="none", lags=0))

# We can still estimate the ECM! but it will just tell us that both
# series are NOT in equilibrium and do not adjust over time

# Engle-Granger
reg <- lm(y~x)
resid <- reg$resid
summary(lm(diff(y) ~ -1 + diff(x) + resid[-t]))


# DHSY
summary(lm(diff(y) ~ diff(x) + y[-t] + x[-t]))


# We can also estimate ECM with stationary data. 
# See e.g. De Boef and Keele, 2008, AJPS







