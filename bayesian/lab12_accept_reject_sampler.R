########################################################################
### Accept-reject sampling                                           ###
### Author: Pablo Barberá                                            ###
### Inspired by:                                                     ###
### http://playingwithr.blogspot.com/2011/06/rejection-sampling.html ###
########################################################################

# Imagine we want to know what the mean and sd of a beta distribution
# with shape parameters a=6 and b=3 are
a <- 6
b <- 3

# This is what the density looks like
x <- seq(0, 1, 0.01)
plot(x, dbeta(x, a, b), type="l", xlab="")

# How would we do it? Luckily, beta distribution has a known analytical
# solution for mean and sd
cat("Mean is equal to", a / (a+b))
cat("Standard deviation is equal to", 
    sqrt((a*b)/((a+b)^2 * (a+b+1))))

# But imagine that we're lazy and do not want to look up formulas on Wikipedia
# Or (as it is often the case), there is no easy analytical solution for these
# moments.

# How would we do it? Monte carlo simulation is one option.

# For example, using ACCEPT REJECT SAMPLING
# with notation from Jackman, p.160, this is the algorithm:

T <- 100000 # number of random draws
theta <- rep(NA, T) # vector to be filled with random draws
c <- 3 # rescaling coefficient

for (t in 1:T){
    z <- runif(1, 0, 1) # majorizing function (uniform from 0 to 1)
    u <- runif(1, 0, 1) # will be used later to accept with probability r
    r <- dbeta(z, a, b) / ( c * dunif(z) ) # accept probability
    if (u <= r){
        theta[t] <- z
    }
    else {
        next
    }
}

# Plotting the results along with the true distribution
hist(theta, breaks=100, freq=FALSE, main="Random draws from beta(6, 3)")
lines(x, dbeta(x, a, b))

# and let's see if we get the same results
cat("Estimated mean is equal to", mean(theta, na.rm=TRUE))
cat("Standard deviation is equal to", sd(theta, na.rm=TRUE))

# Same as above, but more efficient and saving accepted and rejected data to 
# visually see how sampler works

z <- runif(T, 0, 1)
u <- runif(T, 0, 1)
r <- dbeta(z, a, b) / (c * dunif(z) )
accept <- u <= r

d <- data.frame(theta = z, accept = factor(accept, levels=c("TRUE", "FALSE")))

library(ggplot2)
ggplot(d, aes(x=theta, fill=accept)) + geom_histogram(binwidth=0.01)


########################################################################
# Accept-reject sampling from a truncated normal (Jackman page 161)
########################################################################

# density of target function: normal density truncated to region theta > k
p.density <- function(theta, k){
    ifelse(theta>k, exp(-theta^2 / 2), 0)
}

# majorizing function: translated exponential density
m.density <- function(theta, k, alpha){
    alpha * exp( -alpha * (theta - k))
}

# sampling from a standard normal density truncated to lie above k=1

T <- 100000 # number of random draws
theta <- rep(NA, T) # vector to be filled with random draws
alpha <- 1.62 # rescaling coefficient (see p.163)
k <- 1 # truncation point

for (t in 1:T){
    z <- rexp(1, rate=alpha) + k # majorizing function (translated exponential)
    u <- runif(1, 0, 1) # will be used later to accept with probability r
    r <- p.density(z, k) / m.density(z, k, alpha) # accept probability
    if (u <= r){
        theta[t] <- z
    }
    else {
        next
    }
}

# Histogram of sampled values (Figure 3.12 in Jackman)
hist(theta, breaks=100, freq=FALSE, xlim=c(1, 5),
    main="Random draws from standard normal, truncated at 1")
x <- seq(1, 6, 0.01)
lines(x, dnorm(x)/(1-pnorm(1)))

# same as above, but comparing accepted and rejected samples
z <- rexp(T, rate=alpha) + k
u <- runif(T, 0, 1)
r <- p.density(z, k) / m.density(z, k, alpha)
accept <- u <= r
d <- data.frame(theta = z, accept = factor(accept, levels=c("TRUE", "FALSE")))
ggplot(d, aes(x=theta, fill=accept)) + geom_histogram(binwidth=0.01) +
    scale_x_continuous(limits=c(1, 5))





