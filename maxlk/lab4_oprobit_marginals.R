################################################################
## Ordinal probit and marginal effects
## Author: Pablo Barberá
## Quant III Lab 4
## September 24th
################################################################

# installing packages for today
#install.packages("MASS")

# reading the data
data <- read.csv("lab4_data.csv", stringsAsFactors=F)

# Code to run a four-category ordered probit
# Chris Adolph (University of Washington)

# Likelihood for 4 category ordered probit
llk.oprobit4 <- function(param, x, y) {
  # preliminaries
  os <- rep(1, nrow(x))
  x <- cbind(os, x)  
  b <- param[1:ncol(x)]
  t2 <- param[(ncol(x)+1)]
  t3 <- param[(ncol(x)+2)]
  
  # probabilities and penalty function
  xb <- x%*%b
  p1 <- log(pnorm(-xb))
  if (t2<=0)  p2 <- -(abs(t2)*10000)    # penalty function to keep t2>0
  else p2 <- log(pnorm(t2-xb)-pnorm(-xb))
  if (t3<=t2) p3 <- -((t2-t3)*10000)    # penalty to keep t3>t2
  else p3 <- log(pnorm(t3-xb)-pnorm(t2-xb))     
  p4 <- log(1-pnorm(t3-xb)) 

  # -1 * log likelihood (optim is a minimizer)
  -sum(cbind(y==1,y==2,y==3,y==4) * cbind(p1,p2,p3,p4))
}

# Preparing data

y <- data$nytimes   # How often respondent reads NY Times
                    ## 4 = regularly
                    ## 3 = sometimes
                    ## 2 = hardly ever
                    ## 1 = never

x <- cbind(data$age, data$female, data$high.school, data$college)

# Use optim directly
ls.result <- lm(y~x)                    # use ls estimates as starting values
stval <- c(ls.result$coefficients,2,3)  # initial guesses
oprobit.result <- optim(stval, llk.oprobit4, method="BFGS", x=x, y=y, hessian=T)

pe <- oprobit.result$par                # point estimates
vc <- solve(oprobit.result$hessian)     # var-cov matrix
se <- sqrt(diag(vc))                    # standard errors

results <- cbind(pe, se)
dimnames(results) <- list(
    c("Intercept", "age", "female", "high.school", "college", 
        "tau2", "tau3"),
    c("Value", "Std. Error"))
round(results, 5)

# Use MASS polr to do ordered probit
library(MASS)
oprobit <- polr(
    factor(nytimes) ~ age + female + high.school + college, 
    data=data, method="probit")
summary(oprobit)

# predicted probabilities
newdata <- data.frame(age = 18:80, female=mean(data$female),
    high.school=mean(data$high.school), college=mean(data$college))
pred <- predict(oprobit, newdata, type="probs")
str(pred)
plot(18:80, pred[,3])

# using simulation to get standard errors
coefs <- mvrnorm(n=1000, mu=pe, Sigma=vc)
# doing this with a loop, but probably more efficient with *apply
sim.results <- list()
for (age in 18:80){
    x <- c(1, age, mean(data$female), 
        mean(data$high.school), mean(data$college))
    xb <- x %*% t(coefs[,1:5])
    # recall that pr(y=3)=Phi(tau3-xb)-Phi(tau2-xb)
    tau2 <- coefs[,6]; tau3 <- coefs[,7]
    probs <- pnorm(tau3-xb) - pnorm(tau2-xb)
    sim.results[[age]] <- probs
}

sim.results <- do.call(rbind, sim.results)

# getting quantiles
plot.data <- apply(sim.results, 1, function(x) c(quantile(x, .025), mean(x), quantile(x, .975)))

plot(18:80, plot.data[2,], type="l", ylim=c(.05, .25),
    xlab="age", ylab="Pr(y=3)")
lines(18:80, plot.data[1,], lty=3)
lines(18:80, plot.data[3,], lty=3)









