library(gnm)
set.seed(1)
data(occupationalStatus)

RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               Nonlin(MultHomog(origin, destination)), family = poisson,
               data = occupationalStatus)

print(RChomog$deviance, digits=10)
print(RChomog$df)
