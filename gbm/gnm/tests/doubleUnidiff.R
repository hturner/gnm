library(gnm)
set.seed(1)
data(cautres)

doubleUnidiff <- gnm(Freq ~ election*vote + election*class*religion +
                     Mult(Exp(-1 + election), religion:vote) +
                     Mult(Exp(-1 + election), class:vote),
                     family = poisson, data = cautres)

print(doubleUnidiff$deviance, digits=10)
print(doubleUnidiff$df)
