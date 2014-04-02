################################################################
## Probit regression and quantities of interest
## Author: Pablo Barberá
## Quant III Lab 3
## September 19th
################################################################

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
    X <- cbind(1, X)
    loglk <- sum(y * pnorm(X%*%beta, log.p=TRUE) + (1-y) * pnorm(-(X%*%beta), log.p=TRUE))
    return(-loglk)
}

# note: log(1-Phi(X%*%beta)) = log(Phi(-(X%*%beta)))
# why did I do this?

# data
y <- data$twitter
X <- as.matrix(data[,2:6])

# initial values
k <- dim(X)[[2]]+1 ## number of variables: covariates + 1 (constant)
inits <- rnorm(k, 0, 1)

# ML estimation
ml <- optim(par = inits, probit.loglk, X=X, y=y, method="BFGS", hessian = T)
ml$par

# checking with R's built-in function
glm(twitter ~ age + female + high.school + college + nytimes,
    data=data, family=binomial(link="probit"))$coefficients

# computing a few quantities of interest (next recitation we will see more)
probit <- glm(twitter ~ age + female + high.school + college + nytimes,
    data=data, family=binomial(link="probit"))

# predicted probabilities for the average individual
newdata <- data.frame(age = mean(data$age), female=median(data$female),
    high.school=median(data$high.school), college=median(data$college),
    nytimes = median(data$nytimes))
pred <- predict(probit, newdata, type="response", se.fit=TRUE)
pred$fit

# confidence interval
pred$fit - (pred$se.fit * 1.96)
pred$fit + (pred$se.fit * 1.96)

# predicted probabilities at different values of the age variable
newdata <- data.frame(age = 18:80, female=median(data$female),
    high.school=median(data$high.school), college=median(data$college),
    nytimes = median(data$nytimes))
pred <- predict(probit, newdata, type="response", se.fit=TRUE)

plot(18:80, pred$fit, type="l", xlab="age", ylab="Pr(Twitter=1)", ylim=c(0, .40))
lines(18:80, pred$fit - pred$se.fit * 1.96, col="grey80")
lines(18:80, pred$fit + pred$se.fit * 1.96, col="grey80")

# first difference in prob(twitter) by values of nytimes
newdata1 <- data.frame(age = mean(data$age), female=median(data$female),
    high.school=median(data$high.school), college=median(data$college),
    nytimes = 0)
newdata2 <- data.frame(age = mean(data$age), female=median(data$female),
    high.school=median(data$high.school), college=median(data$college),
    nytimes = 1)
pred1 <- predict(probit, newdata2, type="response")
pred2 <- predict(probit, newdata1, type="response")

pred1 - pred2

# note that this is different from the marginal effect of reading the nytimes!
# why? we'll discuss it next time













