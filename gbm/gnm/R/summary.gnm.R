summary.gnm <- function (object, dispersion = NULL, correlation = FALSE,
                          symbolic.cor = FALSE, ...) 
{
    if (is.null(dispersion)) { 
        if (any(object$family$family == c("poisson", "binomial"))) 
            dispersion <- 1
        else if (object$df.residual > 0) {
            if (any(object$weights == 0)) 
                warning("observations with zero weight ",
                        "not used for calculating dispersion")
            dispersion <- sum(object$weights * object$residuals^2)/
                object$df.residual
        }
        else dispersion <- Inf
    }
    if (!"vcov" %in% names(object)){
        start <- coef(object)
        cov.unscaled <- update(object, vcov = TRUE, start = start,
                               trace = FALSE)$vcov
    }
    else cov.unscaled  <- object$vcov
    cov.scaled <- dispersion * cov.unscaled
    ans <- c(object[c("call", "terms", "family", "deviance", "aic",
                      "df.residual", "iter")],
             list(deviance.resid = residuals(object, type = "deviance"), 
                  coefficients = coef(object), dispersion = dispersion,
                  cov.unscaled = cov.unscaled, cov.scaled = cov.scaled))
    if (!is.null(object$original.call))
        ans$original.call <- object$original.call
    if (correlation & object$rank > 0) {
        dd <- sqrt(diag(cov.unscaled))
        ans$correlation <- cov.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.gnm"
    return(ans)
}
