################################################################
## Estimating parameters of poisson/beta distribution with MLE
## Author: Pablo Barberá
## Quant III Lab 1
## September 5th
################################################################

################################################################
## Recovering parameter of a poisson distribution
################################################################

## simulate data for n=100 and lambda=5
y <- rpois(1000, 5)

summary(y)
plot(table(y))

## poisson distribution can only take non-negative integer values

# loglikelihood function
poisson.lk <- function(lambda, y){
    n <- length(y)
    log.lk <- sum(y) * log(lambda) - n * lambda
    return(-log.lk)
}

## last part of log likelihood is dropped because it does not include lambda
## again, also note that we give it a negative sign!
## 'optim' finds the minimum, not the maximum

ml.estimates <- optim(par=1, poisson.lk, y=y, method="BFGS")
ml.estimates$par

# the pdf of the poisson distribution comes with base R
# so we could have simply done:

poisson.lk <- function(lambda, y){
    log.lk <- sum(dpois(y, lambda, log=TRUE))
    return(-log.lk)
}

ml.estimates <- optim(par=1, poisson.lk, y=y, method="BFGS")
ml.estimates$par

# do initial values matter?
ml.estimates <- optim(par=1000, poisson.lk, y=y, method="BFGS")
ml.estimates$par


# assessing convergence of ML estimator using simulation
n <- seq(5, 5000, 25)
convergence <- function(n){
    y <- rpois(n, 5)
    ml.estimates <- optim(par=1, poisson.lk, y=y, method="BFGS")
    return(ml.estimates$par)
}

results <- sapply(n, convergence)

plot(n, results, type="l", ylim=c(4, 6), 
    xlab="Sample size", ylab="Estimated lambda")
abline(h=5, col="grey80")

plot(n, results-5, type="l", ylim=c(-1, 1), 
    xlab="Sample size", ylab="lambda_hat - lambda")
abline(h=0, col="grey80")

# visualizing log.lk distribution
y <- rpois(100, 5)
lambdas <- seq(1, 10, 0.10)
lks <- -sapply(lambdas, poisson.lk, y)

plot(lambdas, lks, type="l", xlab="lambda", ylab="logL")


################################################################
## Recovering parameters of a beta distribution
################################################################

# logLk function for beta distribution
# (different ways of doing the same thing...)

beta.lk <- function(pars, y){
    a <- pars[1]; b <- pars[2]
    log.lk <- -length(y) * 
                log( gamma(a) * gamma(b) / gamma(a+b) ) +
                (a-1) * sum(log(y)) + 
                (b-1)*sum(log(1-y))
    return(-log.lk)
}


beta.lk <- function(pars, y){
    a <- pars[1]; b <- pars[2]
    log.lk <- -length(y) * 
                log( beta(a, b) ) +
                (a-1) * sum(log(y)) + 
                (b-1)*sum(log(1-y))
    return(-log.lk)
}

beta.lk <- function(pars, y){
    a <- pars[1]; b <- pars[2]
    log.lk <- sum(dbeta(y, a, b, log=TRUE))
    return(-log.lk)
}

beta.lk <- function(pars, y) -sum(dbeta(y, pars[1], pars[2], log=TRUE))


# simulation function

simulation <- function(n, a, b){
    y <- rbeta(n, a, b)
    ml.estimates <- optim(c(1,1), beta.lk, y=y)
    return(ml.estimates$par)
}

# Simulating values from a beta distribution and recovering parameters
# using Maximum Likelihood

set.seed(12345) ## note use of seed to ensure replicability
nsizes <- c(10, 100, 1000)
sapply(nsizes, simulation, a=0.50, b=3)

# Changing values to check to what extent it is possible to recover them
mapply(simulation, n=10, a=c(0.50, 5, 10), b=c(3, 10, 50))


# checking consistency of estimator using monte carlo simulations
results <- lapply(seq(10, 1000, 100), 
    function(x) replicate(n=50, simulation(n=x, a=0.50, b=3)))

## the previous line does the same as:
results <- list()
n.sizes <- seq(10, 1000, 100)
for (i in 1:length(n.sizes)){
    results[[i]] <- replicate(n=50, simulation(n=n.sizes[i], a=0.50, b=3))
}

# plotting results
a.est <- unlist(lapply(results, function(x) x[1,]))
b.est <- unlist(lapply(results, function(x) x[2,]))
n.sizes <- rep(seq(10, 1000, 100), each=50)



par (mfrow = c(1,2))
plot(x=n.sizes+runif(length(n.sizes), -50, 50), y=a.est, 
    col = alpha('black', 0.25), pch=16, cex=0.75,
    xlab="Sample Size", ylab="Estimated alpha")
plot(x=n.sizes+runif(length(n.sizes), -50, 50), y=b.est, 
    col = alpha('black', 0.25), pch=16, cex=0.75,
    xlab="Sample Size", ylab="Estimated beta")


# how does value of alpha or beta affect ability to recover parameters?
results <- lapply(c(0.50, seq(1, 1000, 50)), 
    function(x) replicate(n=50, simulation(n=100, a=x, b=3)))
a.est <- unlist(lapply(results, function(x) x[1,]))
alphas <- rep(c(0.50, seq(1, 1000, 50)), each=50)

plot(x=alphas+runif(length(alphas), -25, 25), y=a.est, 
    col = alpha('black', 0.25), pch=16, cex=0.75,
    xlab="True value of alpha", ylab="Estimated alpha")

# why? we need to look at the curvature of the logLk 
contour.beta <- function(a, b){
    y <- rbeta(100, a, b)
    alphas <- seq(0.10, 5.00, 0.05)
    betas <- seq(0.10, 6, 0.05) 
    lks <- matrix(NA, nrow=length(alphas), ncol=length(betas))

    for (i in 1:length(alphas)){
        for (j in 1:length(betas)){
            lks[i,j] <- -beta.lk(pars=c(alphas[i],betas[j]), y=y)
        }
    }   

    lks <- lks/max(lks)
    return(list( alphas = alphas, betas = betas, lks = lks))
}
# levels of contour plot
contours <- c(seq(.05, .95, .10), .99)

par (mfrow = c(1,3))

results <- contour.beta(a=0.50, b=3)
contour(x=results$alphas, y=results$betas, z=results$lks,
    levels=contours, cex=2, xlab="alpha", ylab="beta")

results <- contour.beta(a=1, b=3)
contour(x=results$alphas, y=results$betas, z=results$lks,
    levels=contours, cex=2, xlab="alpha", ylab="beta" )

results <- contour.beta(a=4, b=3)
contour(x=results$alphas, y=results$betas, z=results$lks,
    levels=contours, cex=2, xlab="alpha", ylab="beta")


# extreme example

y <- rbeta(100, 1000, 3)
alphas <- seq(500, 1500, 100)
betas <- seq(0.10, 6, 0.05) 
lks <- matrix(NA, nrow=length(alphas), ncol=length(betas))

for (i in 1:length(alphas)){
    for (j in 1:length(betas)){
        lks[i,j] <- -beta.lk(pars=c(alphas[i],betas[j]), y=y)
    }
}   

lks <- lks/max(lks)
contour(x=alphas, y=betas, z=lks,
    levels=contours, cex=2, xlab="alpha", ylab="beta")


# standard errors

ml.estimates <- optim(c(1,1), beta.lk, y=y, hessian=TRUE)
se <- sqrt(diag(solve(ml.estimates$hessian)))

