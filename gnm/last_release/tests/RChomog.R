library(gnm)
set.seed(1)
data(occupationalStatus)

RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), family = poisson,
               data = occupationalStatus)

print(RChomog$deviance, digits=10)
print(RChomog$df)

data(friend)

###  Fit an association model with homogeneous row-column effects
set.seed(4)
rc2 <- gnm(Freq ~ r + c + Diag(r,c) + instances(MultHomog(r, c), 2),
           family = poisson, data = friend, start = runif(154))

print(rc2$deviance, digits=10)
print(rc2$df)
