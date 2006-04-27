library(gnm)
set.seed(1)
data(yaish)

unidiff <- gnm(Freq ~ educ*orig + educ*dest +
                     Mult(Exp(-1 + educ), orig:dest), family = poisson,
                     data = yaish, subset = (dest != 7))
unidiff
getContrasts(unidiff, grep("\\.1", names(coef(unidiff))))
