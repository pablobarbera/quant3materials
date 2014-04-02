################################################################
## OLS vs Poisson
## Author: Pablo Barberá
## Quant III Lab 6
## October 10th 2013
################################################################

library(foreign)
## Source of data:
# http://emiguel.econ.berkeley.edu/assets/miguel_research/20/Soccer_Replication.zip
data <- read.dta("soccer_data.dta")

## estimating poisson model
summary(
    poisson <- glm(yellow_card ~ civwar + age, family = "poisson", data = data)
    )

## negative binomial
library(MASS)
summary(
    negbin <- glm.nb(yellow_card ~ civwar + age, data = data)
    )


## OLS
summary(
    ols <- lm(yellow_card ~ civwar + age, data = data)
    )


## estimating marginal effects of age using simulation, assuming civwar=0
## POISSON
betas <- mvrnorm(n=1000, mu=poisson$coefficients, 
    Sigma=summary(poisson)$cov.unscaled)
age <- 18:40
marg <- matrix(NA, nrow=1000, ncol=length(age))
for (i in 1:1000){
    for (a in 1:length(age)){
        marg[i,a] <- betas[i,3] * exp(c(1, 0, a) %*% betas[i,])
    }
}

marg <- apply(marg, 2, function(x)
    c(mean(x), quantile(x, c(.025, .975))))

par(mfrow=c(1,3))
plot(age, marg[1,], type="l", 
    xlab="Age", ylab="Marginal effect of age", ylim=c(0.02, 0.07))
lines(age, marg[2,], col="grey")
lines(age, marg[3,], col="grey")

## NEGATIVE BINOMIAL
betas <- mvrnorm(n=1000, mu=negbin$coefficients, 
    Sigma=summary(negbin)$cov.unscaled)
age <- 18:40
marg <- matrix(NA, nrow=1000, ncol=length(age))
for (i in 1:1000){
    for (a in 1:length(age)){
        marg[i,a] <- betas[i,3] * exp(c(1, 0, a) %*% betas[i,])
    }
}

marg <- apply(marg, 2, function(x)
    c(mean(x), quantile(x, c(.025, .975))))

plot(age, marg[1,], type="l", 
    xlab="Age", ylab="Marginal effect of age", ylim=c(0.02, 0.07))
lines(age, marg[2,], col="grey")
lines(age, marg[3,], col="grey")

## OLS
betas <- mvrnorm(n=1000, mu=negbin$coefficients, 
    Sigma=summary(negbin)$cov.unscaled)
age <- 18:40
marg <- matrix(NA, nrow=1000, ncol=length(age))
for (i in 1:1000){
    for (a in 1:length(age)){
        marg[i,a] <- betas[i,3] ## * exp(c(1, 0, a) %*% betas[i,])
    }
}

marg <- apply(marg, 2, function(x)
    c(mean(x), quantile(x, c(.025, .975))))

plot(age, marg[1,], type="l", 
    xlab="Age", ylab="Marginal effect of age", ylim=c(0.02, 0.07))
lines(age, marg[2,], col="grey")
lines(age, marg[3,], col="grey")








