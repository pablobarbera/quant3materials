################################################################
## Identification in logit regression
## Author: Pablo Barberá
## Quant III Lab 5
## October 3th 2013
################################################################

## What happens when we try to estimate a logit/probit model with intercept
## and a threshold parameter?

## log.lk function adding an ancillary parameter
logit.loglk <- function(pars, X, y){
    a <- pars[1]
    b <- pars[2:length(pars)]
    loglk <- sum(y * plogis(a + X%*%b, log.p=TRUE) + (1-y) * plogis(-(a + X%*%b), log.p=TRUE))
    return(-loglk)
}

# reading the data
data <- read.csv("lab3_data.csv", stringsAsFactors=F)

# data
y <- data$twitter
X <- cbind(1, as.matrix(data[,2:3])) ## adding a column of 1's for intercept
k <- dim(X)[[2]] ## number of variables

# 'correct' values
summary(glm(twitter ~ age + female,
    data=data, family=binomial(link="logit")))

# ML estimation
inits <- rep(0, k+1)
logit <- optim(par = inits, logit.loglk, X=X, y=y, method="BFGS", hessian = T)
logit$par

# what happens when we choose different initial values?
inits <- c(-10, 10, rep(0, k-1))
logit <- optim(par = inits, logit.loglk, X=X, y=y, method="BFGS", hessian = T)
logit$par

# what happens to standard errors?
inits <- rep(0, k+1)
logit <- optim(par = inits, logit.loglk, X=X, y=y, method="BFGS", hessian = T)
logit$par
round(sqrt(diag(solve(logit$hessian))), 3)

inits <- c(-10, 10, rep(0, k-1))
logit <- optim(par = inits, logit.loglk, X=X, y=y, method="BFGS", hessian = T)
logit$par
round(sqrt(diag(solve(logit$hessian))), 3)

# another way to understand this problem:
intercepts <- seq(-10, 10, .10)
thresholds <- seq(-10, 10, .10)
lks <- matrix(NA, nrow=length(intercepts), ncol=length(thresholds))

for (i in 1:length(intercepts)){
    for (j in 1:length(thresholds)){
        lks[i,j] <- -logit.loglk(pars=c(intercepts[i], thresholds[j], logit$par[3:4]), X, y)
    }
}   

lks <- max(lks)/lks

contours <- c(seq(.05, .95, .10), .99)
contour(x=intercepts, y=thresholds, z=lks,
    levels=contours, cex=2, xlab="intercept", ylab="threshold")








