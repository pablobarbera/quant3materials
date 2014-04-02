################################################################
## Bayesian Statistics: Introduction to IRT
## Quant III Lab 13
## December 5th 2013
################################################################

library(pscl)

# loading roll call data
# Source: http://jackman.stanford.edu/blog/
load("lab13_senate_rc.rda")
rc <- dropRollCall(rc, dropList=list(codes = "notInLegis", lop = 0))


################################################################
## FITTING IRT MODEL
################################################################

irt <- ideal(rc, store.item=TRUE)

# analysis of results
summary(irt)

# visualizing distribution of ideal points by party
ideal.points <- c(irt$xbar)
party <- rc$legis.data$party
df <- data.frame(ideal.points, party)
library(ggplot2)
ggplot(df, aes(x=ideal.points, fill=party)) + geom_density() + theme_bw() +
    scale_fill_manual(values=c("blue", "red"))

################################################################
## INTERPRETING ITEM PARAMETERS
################################################################


## DIFFICULTY
difficulty <- irt$betabar[,"Difficulty"]
yeas <- as.numeric(rc$vote.data$yeas)
plot(difficulty, yeas)

## DISCRIMINATION
discrimination <- irt$betabar[,"Discrimination D1"]

# top 3 most "discriminatory" bills for POSITIVE values of scale
rc$vote.data[order(discrimination, decreasing=TRUE)[1],]
rc$vote.data[order(discrimination, decreasing=TRUE)[2],]
rc$vote.data[order(discrimination, decreasing=TRUE)[3],]

# top 3 most "discriminatory" bills for NEGATIVE values of scale
rc$vote.data[order(discrimination)[1],]
rc$vote.data[order(discrimination)[2],]
rc$vote.data[order(discrimination)[3],]

# top 5 most informative about ideology
rc$vote.data[order(abs(discrimination), decreasing=TRUE)[1],]
rc$vote.data[order(abs(discrimination), decreasing=TRUE)[2],]
rc$vote.data[order(abs(discrimination), decreasing=TRUE)[3],]

# top 5 least informative about ideology
rc$vote.data[order(abs(discrimination))[1],]
rc$vote.data[order(abs(discrimination))[2],]
rc$vote.data[order(abs(discrimination))[3],]








