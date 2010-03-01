summary.gnm <- function (object, dispersion = NULL, correlation = FALSE,
                         symbolic.cor = FALSE, include.elim = FALSE, ...)
{
    est.disp <- (!object$family$family %in% c("poisson", "binomial") &&
                 is.null(dispersion) && object$df.residual > 0)
    coefs <- coef(object)
    if (include.elim) coefs <- c(object$elim.coefs, coefs)
    if (object$rank > 0) {
        cov.scaled <- vcov(object, dispersion = dispersion,
                           use.eliminate = include.elim)
        ## non-eliminated par only
        estimable <- checkEstimable(object, ...)
        estimable[is.na(estimable)] <- FALSE
        if (is.matrix(cov.scaled))
            sterr <- sqrt(diag(cov.scaled))
        else
            sterr <- diag(cov.scaled)
        is.na(sterr[!estimable]) <- TRUE
        if (include.elim)
            sterr <- c(sqrt(attr(cov.scaled, "varElim")), sterr)
        tvalue <- coefs/sterr
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coefs, sterr, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coefs),
                                         c(dn, "z value", "Pr(>|z|)"))
        }
        else if (object$df.residual > 0) {
            pvalue <- 2 * pt(-abs(tvalue), object$df.residual)
            coef.table <- cbind(coefs, sterr, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coefs),
                                         c(dn, "t value", "Pr(>|t|)"))
        }
        else {
            coef.table <- cbind(coefs, Inf)
            dimnames(coef.table) <- list(names(coefs), dn)
        }
    }
    else {
        coef.table <- matrix(, 0, 4)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
            "t value", "Pr(>|t|)"))
        cov.scaled <- matrix(, 0, 0)
    }
    df.f <- nrow(coef.table)
    non.elim <- seq(object$coef) + nlevels(object$eliminate)*include.elim
    elim <- seq(length.out = length(object$eliminate)*include.elim)
    ans <- c(object[c("call", "ofInterest", "family", "deviance", "aic",
                      "df.residual", "iter")],
             list(deviance.resid = residuals(object, type = "deviance"),
                  coefficients = coef.table[non.elim,],
                  elim.coefs = coef.table[elim,],
                  dispersion = attr(cov.scaled, "dispersion"),
                  df = c(object$rank, object$df.residual, df.f),
                  cov.scaled = as.matrix(cov.scaled)))
    if (correlation & object$rank > 0) {
        dd <- sqrt(diag(cov.scaled))
        ans$correlation <- cov.scaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.gnm"
    ans
}
