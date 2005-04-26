set.seed(1)
data(yaish)

yaish.unidiff <- gnm(Freq ~ educ:orig + educ:dest +
                     Mult(Exp(-1 + educ), orig:dest), family = poisson,
                     data = yaish)
yaish.unidiff
getContrasts(yaish.unidiff, grep("Mult.*educ", names(coef(yaish.unidiff))))
