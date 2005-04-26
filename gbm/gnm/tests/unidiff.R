set.seed(1)
data(yaish)

yaishUnidiff <- gnm(Freq ~ educ:orig + educ:dest +
                     Mult(Exp(-1 + educ), orig:dest), family = poisson,
                     data = yaish)
yaishUnidiff
getContrasts(yaishUnidiff, grep("Mult.*educ", names(coef(yaishUnidiff))))
