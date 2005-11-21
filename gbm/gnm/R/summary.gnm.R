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
        constrain <- object$constrain
        eliminate <- object$eliminate
        needToElim <- seq(sum(!constrain[seq(eliminate)]))[eliminate > 0]
        X <- model.matrix(object)[, !constrain, drop = FALSE]
        Info <- crossprod(X, as.vector(object$weights) * X)
        if (eliminate)
            constrain <- constrain[-seq(eliminate)]
        if (sum(constrain) > 0) {
            cov.unscaled <- array(0, dim = rep(length(theta), 2),
                              dimnames = rep(list(names(theta)), 2))
            cov.unscaled[!constrain, !constrain] <-
                MPinv(Info, eliminate = needToElim, onlyNonElim = TRUE)
        }
        else
            cov.unscaled <- MPinv(Info, eliminate = needToElim,
                                  onlyNonElim = TRUE)
        attr(cov.unscaled, "rank") <- NULL
    }
    else cov.unscaled  <- object$vcov
    cov.scaled <- dispersion * cov.unscaled
    ans <- c(object[c("call", "terms", "family", "deviance", "aic",
                      "df.residual", "iter")],
             list(deviance.resid = residuals(object, type = "deviance"), 
                  coefficients = coef(object), dispersion = dispersion,
                  cov.unscaled = cov.unscaled, cov.scaled = cov.scaled))
    if (correlation & object$rank > 0) {
        dd <- sqrt(diag(cov.unscaled))
        ans$correlation <- cov.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.gnm"
    return(ans)
}
