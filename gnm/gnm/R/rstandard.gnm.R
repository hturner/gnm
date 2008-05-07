rstandard.gnm <- function(model, ...) {
    so <- summary(model)
    res <- na.omit(so$deviance.resid[model$prior.weights != 0])
    res <- naresid(model$na.action, res)
    res <- res/sqrt(so$dispersion * (1 - hatvalues(model)))
    res[is.infinite(res)] <- NaN
    if (!is.null(model$table.attr))
        attributes(res) <- model$table.attr
    res
}
