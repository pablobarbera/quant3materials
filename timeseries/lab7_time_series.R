################################################################
## Intro to time-series with R
## Author: Pablo Barberá
## Quant III Lab 7
## October 17th 2013
################################################################

library(car)
set.seed(123)

################################################################
## TIME-SERIES WITH IID ERRORS AND NO DISTRIBUTED LAGS
## (i.e., no serial correlation and no lagged covariates)
################################################################

# generating simulated data
t <- 200
x <- runif(t, -1, 1)
e <- rnorm(t, 0, 1)
beta <- 3

y <- beta * x + e

plot(ts(y))


# OLS is unbiased and efficient
summary(ols <- lm(y ~ x))

# testing for serial correlation
durbinWatsonTest(ols, max.lag=1)
# note: H0 = no autocorrelation
# D-W Stats near 2 imply we do not reject the null

# LM test = Breusch-Godfrey
summary(ols <- lm(ols$resid ~ x + c(NA, ols$resid[-t])))

Rstat <- (t - 1) * summary(ols)$r.squared
Rstat
qchisq(0.95, 1) ## critical value
p <- pchisq(Rstat, df=1) ## p-value: cannot reject null no serial correlation
1-p ## p>0.05, so we do not reject the null of no serial correlation

# IMPULSE RESPONSE FUNCTION for x
# (x goes up one unit in t and then back down to zero in t+1)
new.x <- x
new.x[30] <- new.x[30] + 1 ## ading one unit to obs. 300
new.y <- beta * new.x + e

plot(1:t, y, type="l", xlim=c(25, 35))
lines(1:t, new.y, col="blue")

# easier to see without error and assuming covariate is at some eq. value
# (see Beck and Katz, 2011, for more examples)
x <- rep(5, t)
new.x <- x; new.x[30] <- new.x[30] + 1

y <- beta * x
new.y <- beta * new.x 

plot(1:t, y, type="l", xlim=c(25, 35))
lines(1:t, new.y, col="blue")


# UNIT RESPONSE FUNCTION for x
# (x goes up one unit in t and remains up one unit for all remaining periods)
new.x <- x; new.x[30:t] <- new.x[30:t] + 1
new.y <- beta * new.x 

plot(1:t, y, type="l", xlim=c(25, 35))
lines(1:t, new.y, col="blue")


################################################################
## TIME-SERIES WITH IID ERRORS AND DISTRIBUTED LAGS
## (i.e., no serial correlation and lagged covariates)
################################################################

rm(list=ls())
set.seed(777)

t <- 200
## now using 'ts' function (this will make computing lags easier)
x <- ts(runif(t, -1, 1))
l.x <- lag(x, -1)
e <- ts(rnorm(t, 0, 1))
beta <- 3
gamma <- 1
y <- beta * x + gamma * l.x + e

## putting it all together
d <- ts.union(y, x, l.x, e)
plot(d)

# everything still works...
summary(ols <- lm(y ~ x + l.x, data=d))

durbinWatsonTest(ols, max.lag=1) ## no serial correlation

# same with LM test
resid <- ts(ols$resid, start=2)
l.resid <- lag(resid, -1)
d <- ts.union(y, x, l.x, resid, l.resid)
summary(ols <- lm(resid ~ x + l.x + l.resid, data=d))
Rstat <- (t - 1) * summary(ols)$r.squared
qchisq(0.95, 1) ## critical value
p <- pchisq(Rstat, df=1) ## p-value: cannot reject null no serial correlation
1-p ## p>0.05, so we do not reject the null of no serial correlation

# IMPULSE RESPONSE FUNCTION for x
# (x goes up one unit in t and then back down to zero in t+1)
x <- ts(rep(5, t)); new.x <- x; new.x[30] <- new.x[30] + 1
l.x <- lag(x, -1); new.l.x <- lag(new.x, -1)
y <- beta * x + gamma * l.x
new.y <- beta * new.x + gamma * new.l.x

plot(y, xlim=c(25, 35))
lines(new.y, col="blue")
diff <- new.y - y
diff[25:35] ## in t, it goes up by beta; in t+1, it goes up by gamma

# UNIT RESPONSE FUNCTION for x
# (x goes up one unit in t and remains up one unit for all remaining periods)
x <- ts(rep(5, t)); new.x <- x; new.x[30:t] <- new.x[30:t] + 1
l.x <- lag(x, -1); new.l.x <- lag(new.x, -1)
y <- beta * x + gamma * l.x
new.y <- beta * new.x + gamma * new.l.x

plot(y, xlim=c(25, 35))
lines(new.y, col="blue")
diff <- new.y - y
diff[25:35] ## in t, it goes up by beta; in t+1..., it goes up by beta+gamma


################################################################
## TIME-SERIES WITH AR1 ERRORS
## (i.e., serial correlation and no lagged covariates)
################################################################

rm(list=ls())
set.seed(444)

t <- 200
x <- ts(runif(t, -1, 1))
rho <- 0.2
# generating AR1 errors
nu <- rnorm(t)
e <- c(nu[1], rep(NA, t-1))
for (i in 2:t){
    e[i] <- rho * e[i-1] + nu[i]
}
e <- ts(e)
beta <- 3
y <- beta * x + e


## putting it all together
d <- ts.union(y, x, e)
plot(d) ## what plots show serial correlation?

# unbiased BUT inefficient
summary(ols <- lm(y ~ x, data=d))

durbinWatsonTest(ols, max.lag=1) ## we reject the null of no serial corr.

resid <- ts(ols$resid, start=2)
l.resid <- lag(resid, -1)
d <- ts.union(y, x, resid, l.resid)

summary(ols <- lm(resid ~ x + l.resid, data=d))
Rstat <- (t - 1) * summary(ols)$r.squared
qchisq(0.95, 1) ## critical value
p <- pchisq(Rstat, df=1) ## p-value: we reject null no serial correlation
1-p ## p<0.05, so now we DO reject the null of no serial correlation

## how to solve this?
## One option is to fix standard errors using prais-winsten or cochrane-orcutt
## Another option is to add lagged DVs, which (usually) solves the problem:

l.y <- lag(y, -1)
d <- ts.union(y, x, l.y, e)
summary(ols <- lm(y ~ l.y + x, data=d))

## (remember that DW does not work when you have lagged DV)
## (also note that if you do NOT have serial correlation and use a lagged DV,
## then estimates are biased! See e.g. Keele and Kelly, 2006)

resid <- ts(ols$resid, start=2) ## 
d <- ts.union(y, x, resid, l.resid, l.y)

summary(ols <- lm(resid ~ x + l.resid + l.y, data=d))
Rstat <- (t - 1) * summary(ols)$r.squared
qchisq(0.95, 1) ## critical value
p <- pchisq(Rstat, df=1) ## p-value: still reject null no serial correlation
1-p ## p<0.05, so now we DO reject the null of no serial correlation


# IMPULSE RESPONSE FUNCTION for nu (iid error)
# (x goes up one unit in t and then back down to zero in t+1)
## 1) baseline assuming x and nu are at equilibrium values
x <- ts(rep(5, t))
nu <- rep(0, t) 
e[1] <- 0
for (i in 2:t){
    e[i] <- rho * e[i-1] + nu[i]
}
e <- ts(e)
y <- beta * x + e

## 2) new values, after adding +1 to nu for one period
nu[30] <- nu[30] + 1
for (i in 2:t){
    e[i] <- rho * e[i-1] + nu[i]
}
e <- ts(e)
new.y <- beta * x + e

plot(y, xlim=c(25, 35))
lines(new.y, col="blue")
diff <- new.y - y
diff[25:35] 

################################################################
## what happens if rho is now .9 ?
################################################################

rho <- .9

x <- ts(rep(5, t))
nu <- rep(0, t) 
e[1] <- 0
for (i in 2:t){
    e[i] <- rho * e[i-1] + nu[i]
}
e <- ts(e)
y <- beta * x + e

nu[30] <- nu[30] + 1
for (i in 2:t){
    e[i] <- rho * e[i-1] + nu[i]
}
e <- ts(e)
new.y <- beta * x + e

plot(y, xlim=c(25, 50))
lines(new.y, col="blue")
diff <- new.y - y
diff[25:50] 

################################################################
## what if rho = 1? (random walk)
################################################################

rm(list=ls())
set.seed(333)

t <- 200
x <- ts(runif(t, -1, 1))
rho <- 1
# generating AR1 errors
nu <- rnorm(t)
e <- c(nu[1], rep(NA, t-1))
for (i in 2:t){
    e[i] <- rho * e[i-1] + nu[i]
}
e <- ts(e)
beta <- 3
y <- beta * x + e


## putting it all together
d <- ts.union(y, x, e)
plot(d)  

# everything is still fine
summary(ols <- lm(y ~ x, data=d))



################################################################
## what if rho > 1? (non-stationary time-series)
################################################################

rm(list=ls())
set.seed(333)

t <- 200
x <- ts(runif(t, -1, 1))
rho <- 1.05
# generating AR1 errors
nu <- rnorm(t)
e <- c(nu[1], rep(NA, t-1))
for (i in 2:t){
    e[i] <- rho * e[i-1] + nu[i]
}
e <- ts(e)
beta <- 3
y <- beta * x + e


## putting it all together
d <- ts.union(y, x, e)
plot(d)  

# everything explodes
summary(ols <- lm(y ~ x, data=d))









