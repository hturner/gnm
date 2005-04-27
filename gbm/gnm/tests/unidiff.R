library(gnm)
set.seed(1)
data(yaish)

unidiff <- gnm(Freq ~ educ:orig + educ:dest +
                     Mult(Exp(-1 + educ), orig:dest), family = poisson,
                     data = yaish)
unidiff
getContrasts(unidiff, grep("Mult.*educ", names(coef(unidiff))))
