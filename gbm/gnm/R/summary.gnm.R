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
    if (!("vcov" %in% names(object))){
        eliminate <- object$eliminate
        na <- is.na(coef(object)) ## either constrained or aliased
        needToElim <- seq(sum(!na[seq(eliminate)]))[eliminate > 0]
        X <- model.matrix(object)[, !na, drop = FALSE]
        Info <- crossprod(X, as.vector(object$weights) * X)
        if (sum(na) > 0) {
            cov.unscaled <- array(0, dim = rep(length(na), 2),
                              dimnames = rep(list(names(coef(object))), 2))
            cov.unscaled[!na, !na] <-
                MPinv(Info, eliminate = needToElim, onlyNonElim = FALSE)
        }
        else
            cov.unscaled <- MPinv(Info, eliminate = needToElim,
                                  onlyNonElim = FALSE)
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
