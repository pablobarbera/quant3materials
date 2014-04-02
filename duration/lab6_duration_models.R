################################################################
## Duration models with R
## Author: Pablo Barberá
## Quant III Lab 6
## October 10th 2013
################################################################

library(survival)
library(KMsurv)

data <- read.csv("billboard_data.csv", stringsAsFactors=F)
data <- na.omit(data)

################################################################
## NON-PARAMETRIC ANALYSIS
################################################################

# density of failure times, f(t)
plot(density(data$duration))

# failure function F(t) or CDF for f(t)
fit <- survfit(Surv(duration, event=failure) ~ 1, data=data)
plot(fit, xlab="Weeks", ylab="F(t)", fun="event")
# (note that first we need to fit a model with just a constant)

# survivor function S(t) = 1 - F(t) OR CDF for (1-f(t))
plot(fit, xlab="Weeks", ylab="S(t)")

# hazard rate h(t) = f(t) / S(t)  OR  f(t) / (1 - F(t))
# discrete hazard rate (note that Stata applies some smoothing to this)
library(muhaz)
fit2 <- kphaz.fit(data$duration, data$failure)
kphaz.plot(fit2)


#############     PARAMETRIC MODELS      #######################

################################################################
# EXPONENTIAL DISTRIBUTION: CONSTANT HAZARD RATE
################################################################

# Estimation

exponential <- survreg(
    Surv(duration, event=failure) ~ 1 + year + numberone + length_song + bpm,
    data=data, dist="exponential")

summary(exponential)
round(exponential$coefficients, 3)

# Marginal effects
x <- c(1, 1995, 0, 3, mean(data$bpm))
b <- exponential$coefficients
b[4] * exp(x %*% b)

# (and you can use simulation to get standard errors)

# baseline hazard rate
t <- seq(0, 70, 1)
lambda.i <- exp(-predict(exponential, data, type="linear"))
lambda <- mean(lambda.i, na.rm=TRUE)

hazard <- lambda
plot(t, rep(hazard, length(t)), type="l", main="Exponential", xlab="Weeks in Billboard Chart",
   ylab="Hazard Rate")


################################################################
# WEIBULL DISTRIBUTION: MONOTONICALLY INCREASING/DECREASING HAZARD RATE
################################################################

weibull <- survreg(
    Surv(duration, event=failure) ~ 1 + year + numberone + length_song + bpm,
    data=data, dist="weibull")

summary(weibull)
round(weibull$coefficients, 3)

# Marginal effects
x <- c(1, 1995, 0, 3, mean(data$bpm))
b <- weibull$coefficients
b[4] * exp(x %*% b)

# Hazard rate

lambda.i <- exp(-predict(weibull, data, type="linear"))
lambda <- mean(lambda.i, na.rm=TRUE)
t <- seq(0,70,1)
p <- 1/weibull$scale
scale <- weibull$scale
hazard <- lambda * p * (lambda * t)^(p-1)
plot(t, hazard, type="l", main="Weibull", 
    xlab="Weeks in Billboard Chart", ylab="Hazard Rate")


################################################################
# LOGNORMAL DISTRIBUTION: NONMONOTONIC  HAZARD RATE
################################################################

lognormal <- survreg(
    Surv(duration, event=failure) ~ 1 + year + numberone + length_song + bpm,
    data=data, dist="lognormal")

summary(lognormal)
round(lognormal$coefficients, 3)

# Marginal effects
x <- c(1, 1995, 0, 3, mean(data$bpm))
b <- lognormal$coefficients
b[4] * exp(x %*% b)

# Hazard Rate in Lognormal
lambda.i <- exp(-predict(lognormal, data, type="linear"))
lambda <- mean(lambda.i, na.rm=TRUE)
p <- 1/exponential$scale
pdf <- (2*pi)^{-1/2} * p * t^{-1} * exp((-p^2 * (log(lambda*t))^2)/2)
cdf <- 1 - pnorm(p*log(lambda*t))

hazard <- pdf/cdf
plot(t, hazard, type="l", main="Log-normal", 
    xlab="Weeks in Billboard Chart", ylab="Hazard Rate")


################################################################
# COX MODEL
################################################################

data <- read.dta("civil_cox.dta")
names(data)[36:37] <- c("t", "t0")

(cox.model <- coxph(Surv(date0, date1, event=cens, type="counting") ~ 
    gini_m + ginmis + rgdpch + elf + elf2 + logpop + y70stv + y80stv + 
    y90stv + d2 + d3 + d4 + cluster(indsp), data=data, method="efron"))

hazard <- basehaz(cox.model)
plot(hazard$time, hazard$hazard, type="l", ylab="baseline hazard", xlab="time")



