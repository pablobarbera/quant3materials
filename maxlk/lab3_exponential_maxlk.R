################################################################
## Estimating lambda of an exponential distribution
## Author: Pablo Barberá
## Quant III Lab 3
## September 19th
################################################################

################################################################
## Plotting the likelihood function
################################################################

# As usual, different ways of doing this

exp.lk <- function(lambda, y){
    lk <- prod(lambda * exp(-lambda * y))
    return(lk)
}

exp.lk <- function(lambda, y){
    lk <-prod(dexp(y, lambda))
    return(lk)
}

exp.lk <- function(lambda, y) prod(dexp(y, lambda))

# now we generate simulated data from the exponential
# (note: lambda=3 and n=100)
y <- rexp(100, 3)

# and then we plot the likelihood
# (true value of lambda is 3, and then we explore a 'grid' of possible lambdas)
lambdas <- seq(0, 5, 0.10)
lks <- sapply(lambdas, exp.lk, y)
plot(lambdas, lks, type="l", xlab="lambda", ylab="likelihood")

# note that this is the same as:
lks <- c()
for (lambda in lambdas){
    lks <- c(lks, exp.lk(lambda, y))
}
plot(lambdas, lks, type="l", xlab="lambda", ylab="likelihood")


# now we replicate this for different samples sizes
exp.sim <- function(n, true.lambda){
    y <- rexp(n, true.lambda)
    lambdas <- seq(0, 5, 0.10)
    lks <- sapply(lambdas, exp.lk, y)
    # normalizing (not really necessary)
    lks <- lks/max(lks)
    return(
        list(data.frame(
            size=paste0("n=", n),
            lambdas=lambdas,
            lks=lks)))
}

set.seed(777)
sim.values <- sapply(c(10, 50, 100, 500), exp.sim, true.lambda=3)

# note that the function will return a list of data frames
# we can plot them individually...
plot(sim.values[[4]]$lambdas, sim.values[[4]]$lks, 
    type="l", xlab="lambda", ylab="likelihood")


# or do a slightly more sophisticated plot
plot.data <- do.call(rbind, sim.values)
library(ggplot2)
plot <- ggplot(plot.data, aes(x=lambdas, y=lks, group=size)) +
        geom_line(aes(color=size)) +
        scale_color_discrete("Sample Size") +
        scale_x_continuous("lambda") +
        scale_y_continuous("Likelihood") +
        geom_vline(xintercept=3, color="grey70") +
        theme_bw() +
        theme(
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
plot
ggsave("likelihood_exp.pdf", plot, width=6, height=4)


################################################################
## Maximum likelihood for an exponential distribution
################################################################

# suppose we choose true lambda = 2
n <- 100; true.lambda <- 2
y <- rexp(n, true.lambda)
# then mean(y)=1/2
mean(y)
# and when we estimate lambda hat, its standard error should be:
true.lambda / sqrt(n)
# i.e., 0.20
# in this sample, we should get:
lambda.hat <- 1/mean(y)
lambda.hat / sqrt(n)

# we already know how this goes...

loglk.exp <- function(lambda, y){
    -sum(dexp(y, lambda, log=TRUE))
}

simulation <- function(n, true.lambda){
    y <- rexp(n, true.lambda)
    ml <- optim(1, loglk.exp, method="BFGS", y=y, hessian=TRUE)
    return(c(ml$par, ml$hessian))
}

est <- replicate(n=1000, simulation(n=100, true.lambda=true.lambda))

# histogram to see distribution of lambda hat
hist(est[1,], xlab="lambda hat", main=NULL, xlim=c(1, 3))

# and histogram to see distribution of s.e. of lambda hat
hist(sqrt(1/est[2,]), xlab="s.e. lambda hat", main=NULL, xlim=c(0.10, 0.30))

# note what we are doing:
# for each simulation, we compute sqrt(1/hessian):
# in the first simulation:
est[,1]
# where first element is lambda hat, second is hessian
# so s.e. of lambda hat in that specific simulation is:
hessian <- est[2,1]
sqrt( 1 / hessian)


# now we can do it again for sample size=1600, where s.e. lambda hat
# should be:
true.lambda / sqrt(1600)
# and compare

par (mfrow = c(2,2), mar=c(3,3,2,1), mgp=c(2,.7,0), bty="l", tck=-.02)
# sample size=100
est <- replicate(n=1000, simulation(n=100, true.lambda=true.lambda))
hist(est[1,], xlab="lambda hat", main=NULL, xlim=c(1, 3))
hist(sqrt(1/est[2,]), xlab="s.e. lambda hat", main=NULL, xlim=c(0.10, 0.30))
# sample size=1600
est <- replicate(n=1000, simulation(n=1600, true.lambda=true.lambda))
hist(est[1,], xlab="lambda hat", main=NULL, xlim=c(1, 3))
hist(sqrt(1/est[2,]), xlab="s.e. lambda hat", main=NULL, xlim=c(0.0, 0.10))





