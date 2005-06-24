library(gnm)
set.seed(1)
data(voting)

count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)

classMobility <- gnm(yvar ~ Nonlin(Dref(origin, destination)), 
                     family = binomial, data = voting)
classMobility

upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

socialMobility <- gnm(yvar ~ Nonlin(Dref(origin, destination,
                                         formula = ~ 1 + downward + upward),
                                    data = data.frame(downward, upward)),
                      family = binomial, data = voting)
socialMobility
