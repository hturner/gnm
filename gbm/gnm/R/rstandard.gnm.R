rstandard.gnm <- function(object) {
    so <- summary(object)
    res <- na.omit(so$deviance.resid[object$prior.weights != 0])
    res <- naresid(object$na.action, res)
    res/sqrt(so$dispersion * (1 - hatvalues(object)))
}
