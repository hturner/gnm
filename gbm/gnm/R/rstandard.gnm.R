rstandard.gnm <- function(model, ...) {
    so <- summary(model)
    res <- na.omit(so$deviance.resid[model$prior.weights != 0])
    res <- naresid(model$na.action, res)
    res/sqrt(so$dispersion * (1 - hatvalues(model)))
}
