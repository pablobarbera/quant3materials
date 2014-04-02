################################################################
## 4 things you might want to know about R
## Author: Pablo Barberá
## Quant III Lab 1
## September 5th
################################################################

################################################################
## 1) R functions are like math functions
################################################################

## a trivial example
square_prod <- function(x,y){
    return((x * y)^2)
}

square_prod(2, 1)
square_prod(5, 10)

## computing OLS estimates with matrix algebra
ols <- function(y, X){
    if (all(X[1,]!=1)){ X <- cbind(1, X) }
    beta <- solve(t(X)%*%X)%*%t(X)%*%y
    rownames(beta)[1] <- "(Intercept)"
    return(round(t(beta),3))
}

## simulating data
x1 <- rnorm(100, 0, 1)
x2 <- rnorm(100, 0, 1)
y <- 5 + 3 * x1 + -1 * x2 + rnorm(100, 0, 1)
X <- cbind(x1, x2)

## running function and checking with R built-in estimator
ols(y, X)
lm(y ~ x1 + x2)


################################################################
## 2) *apply functions are just loops
################################################################

## 'sapply' applies a function for each element of a vector

random_numbers <- rnorm(10, 0, 1)

# trivial example: add +1 to each element of a vector
sapply(random_numbers, FUN=function(x) x+1 )

# same as:
for (r in random_numbers){ print(r + 1) }

# but note that many functions in R are vectorized
random_numbers + 1

# believe it or not, the '+' sign is also a function in R!
# (try typing ?'+')

## 'apply' works with matrices and applies a function to each
## row or column of that matrix

mat <- matrix(rnorm(100), nrow=10, ncol=10)

# trivial example: compute the mean of the values in each row
apply(mat, 1, mean)

# same but for columns
apply(mat, 2, mean)

# loop equivalent
for (row in 1:dim(mat)[1]){
    print( mean(mat[row,]) )
}

# and the vectorized equivalent
rowMeans(mat)

## there are others: lapply, vapply, mapply...
## and a package, plyr, that adds many more
## (we might deal with them later in the semester)


################################################################
## 3) Use str to learn more about objects
################################################################

str(mat)

data <- read.csv("lab1_data.csv")
str(data)


################################################################
## 4) stringsAsFactors=F should be the default
################################################################

data$tweets[1]
data$tweets <- as.numeric(data$tweets)
data$tweets[1]

# Ugh.

data <- read.csv("lab1_data.csv", stringsAsFactors=F)
data$tweets <- as.numeric(data$tweets)
data$tweets[1]

summary(data$tweets)



