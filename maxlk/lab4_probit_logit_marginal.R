################################################################
## Probit / Logit and marginal effects
## Author: Pablo Barberá
## Quant III Lab 4
## September 24th
################################################################

# installing packages for today
#install.packages("MASS")

# reading the data
data <- read.csv("lab3_data.csv", stringsAsFactors=F)

str(data)
summary(data)
# 2901 observations (with no missing values)
# 6 variables:
# twitter -- whether respondent has an active Twitter accounts
# age
# female
# high.school -- whether respondent finished grad school
# college -- whether respondent has a college degree
# (reference category for these two is no education/primary schoool)
# nytimes: whether respondent reads ny times

# we want to see what variables predict twitter use

# log-likelihood function

probit.loglk <- function(beta, X, y){
    loglk <- sum(y * pnorm(X%*%beta, log.p=TRUE) + (1-y) * pnorm(-(X%*%beta), log.p=TRUE))
    return(-loglk)
}

# note: log(1-Phi(X%*%beta)) = log(Phi(-(X%*%beta)))
# why did I do this?

logit.loglk <- function(beta, X, y){
    loglk <- sum(y * plogis(X%*%beta, log.p=TRUE) + (1-y) * plogis(-(X%*%beta), log.p=TRUE))
    return(-loglk)
}

logit.loglk <- function(beta, X, y){
    loglk <- sum(-log(1+exp(X%*%beta))+y * X%*%beta  )
    return(-loglk)
}

# data
y <- data$twitter
X <- cbind(1, as.matrix(data[,2:6])) ## adding column of 1's for intercept

# initial values
k <- dim(X)[[2]] ## number of variables
inits <- rnorm(k, 0, 1)

# ML estimation
probit <- optim(par = inits, probit.loglk, X=X, y=y, method="BFGS", hessian = T)
probit$par
sqrt(diag(solve(probit$hessian)))

# checking with R's built-in function
glm(twitter ~ age + female + high.school + college + nytimes,
    data=data, family=binomial(link="probit"))$coefficients

# Same for logit
logit <- optim(par = inits, logit.loglk, X=X, y=y, method="BFGS", hessian = T)
logit$par
sqrt(diag(solve(logit$hessian)))

glm(twitter ~ age + female + high.school + college + nytimes,
    data=data, family=binomial(link="logit"))$coefficients

# average marginal effect of college
xb <- X %*% probit$par
marg <- dnorm(xb) * probit$par[5]
mean(marg)

# marginal effect of college for average individual
xb <- colMeans(X) %*% probit$par
marg <- dnorm(xb) * probit$par[5]
marg

# marginal effect of college for female age=25 nytimes=0
xb <- c(1, 25, 0, 0, mean(data$college), 0) %*% probit$par
marg <- dnorm(xb) * probit$par[5]
marg

# predicted probabilities for college=0 and college=1
xb <- c(1, mean(data$age), mean(data$female), 
    mean(data$high.school), 0, mean(data$nytimes)) %*% probit$par
pred <- pnorm(xb)
pred

xb <- c(1, mean(data$age), mean(data$female), 
    mean(data$high.school), 1, mean(data$nytimes)) %*% probit$par
pred <- pnorm(xb)
pred

# using simulation to compute standard errors for marginal effects
library(MASS)
# simulate 1000 draws from coefficient vector
sims <- 1000
coefs <- mvrnorm(n=sims, mu=probit$par, Sigma=solve(probit$hessian))

# inefficient way of doing this:
results <- c()
for (i in 1:sims){
    xb <- X %*% coefs[i,]
    marg <- dnorm(xb) * coefs[i, 5] ## note that uncertainty enters here too
    results <- c(results, mean(marg))
}
# estimate of marginal effect
mean(results)
# should be the same as...
xb <- X %*% probit$par
mean(dnorm(xb) * probit$par[5])
# and now we can get quantiles to get 95% confidence interval
quantile(results, probs=c(.025, .975))

# faster with matrix algebra and *apply
xb <- X %*% t(coefs)
marg <- sapply(1:1000, function(x) dnorm(xb[,x]) * coefs[x,5])
margs <- colMeans(marg)
mean(margs)
quantile(margs, probs=c(.025, .975))

# now the same for predicted values
xb <- cbind(1, 18:80, mean(data$female), mean(data$high.school), 
    mean(data$college), mean(data$nytimes)) %*% t(coefs)
pred <- pnorm(xb)
preds <- rowMeans(pred)
preds.ci <- apply(pred, 1, quantile, prob=c(0.025, 0.975))

# plot showing predicted probabilities, with CI
plot(18:80, preds, type="l", ylim=c(0, .40),
    xlab="age", ylab="Pr(y=1)")
lines(18:80, preds.ci[1,], lty=3)
lines(18:80, preds.ci[2,], lty=3)

# sample plot with ggplot2
plot.data <- data.frame(
    age=18:80, pred = preds, lo = preds.ci[1,], hi=preds.ci[2,])

library(ggplot2)
plot <- ggplot(plot.data, aes(x=age, y=pred)) +
        geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.25) + 
        geom_line() +
        scale_x_continuous(limits=c(20,80)) +
        scale_y_continuous("Predicted Prob. of Twitter Use", limits=c(0, .35)) +
        theme_bw() +
        theme(
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 
plot

# what about first differences? 
# change in pred. prob from college=0 to college=1
xb1 <- cbind(1, mean(data$age), mean(data$female), mean(data$high.school), 
    0, mean(data$nytimes)) %*% t(coefs)
xb2 <- cbind(1, mean(data$age), mean(data$female), mean(data$high.school), 
    1, mean(data$nytimes)) %*% t(coefs)
fd <- pnorm(xb2) - pnorm(xb1)

mean(fd)
quantile(fd, probs=c(0.025, 0.975))





