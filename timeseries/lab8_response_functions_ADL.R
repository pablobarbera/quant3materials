################################################################
## Impulse/Unit response functions for ADL model
## Author: Pablo Barberá
## Quant III Lab 8
## October 24th 2013
################################################################

rm(list=ls())
set.seed(123)

# ADL model:
# y_t = phi * y_t-1 + beta * x + gamma * x_t-1 + epsilon

# parameters
phi <- 0.95
beta <- 3
gamma <- 2
t <- 200

response.function <- function(phi, beta, gamma, t=200, impulse=FALSE, unit=FALSE){
    # preparing covariates
    x <- rep(0, t)
    if (impulse) x[100] <- x[100] + 1
    if (unit) x[100:t] <- x[100:t] + 1
    l.x <- c(0, x[-t])
    # preparing dv
    y <- rep(NA, t)
    y[1] <- phi * 0 + x[1] * beta + l.x[1] * gamma 
    for (i in 2:t){
        y[i] <- phi * y[i-1] + x[i] * beta + l.x[i] * gamma 
    }
    rng <- c(min(y[95:110])-1,  max(y[95:110])+1)
    par(mar=c(3, 3, 2, 1.5))

    plot((95:110)-99, y[95:110], ylim=rng)
    lines((95:110)-99, y[95:110])

    mtext("t", side=1, line=2)
    mtext("y", side=2, line=2)
    if (impulse) mtext("Impulse Response Function", side=3, line=0.2, cex=0.8)
    if (unit) mtext("Unit Response Function", side=3, line=0.2, cex=0.8)

}


# impulse response functions when phi = 0
response.function(phi=0, beta, gamma, impulse=TRUE, plot=4)

# unit response functions when phi = 0
response.function(phi=0, beta, gamma, unit=TRUE, plot=4)


# impulse response functions when phi = 1
response.function(phi=1, beta, gamma, impulse=TRUE, plot=4)

# unit response functions when phi = 1
response.function(phi=1, beta, gamma, unit=TRUE, plot=4)

# impulse response functions when gamma=0
response.function(phi=0.95, beta, gamma=0, impulse=TRUE, plot=4)

# unit response functions when gamma=0
response.function(phi=0.95, beta, gamma=0, unit=TRUE, plot=4)

# impulse response functions when beta=0
response.function(phi=0.95, beta=0, gamma, impulse=TRUE, plot=4)

# unit response functions when beta=0
response.function(phi=0.95, beta=0, gamma, unit=TRUE, plot=4)

# impulse response functions when gamma = - phi * beta
response.function(phi=0.95, beta, gamma=-phi * beta, impulse=TRUE, plot=4)

# unit response functions when gamma = - phi * beta
response.function(phi=0.95, beta, gamma=-phi * beta, unit=TRUE, plot=4)









