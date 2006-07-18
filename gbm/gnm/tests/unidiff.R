library(gnm)
set.seed(1)
data(yaish)

unidiff <- gnm(Freq ~ educ*orig + educ*dest +
                     Mult(Exp(-1 + educ), orig:dest), family = poisson,
                     data = yaish, subset = (dest != 7))

print(unidiff$deviance, digits = 10)
print(unidiff$df)

getContrasts(unidiff, grep("Mult.*educ", names(coef(unidiff))))
