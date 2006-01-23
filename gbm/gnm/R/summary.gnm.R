summary.gnm <- function (object, dispersion = NULL, correlation = FALSE,
                          symbolic.cor = FALSE, ...) 
{
    est.disp <- FALSE
    if (is.null(dispersion)) { 
        if (any(object$family$family == c("poisson", "binomial"))) 
            dispersion <- 1
        else if (object$df.residual > 0) {
            est.disp <- TRUE
            if (any(object$weights == 0)) 
                warning("observations with zero weight ",
                        "not used for calculating dispersion")
            dispersion <- sum(object$weights * object$residuals^2)/
                object$df.residual
        }
        else dispersion <- Inf
    }
    coefs <- coef(object)
    if (object$rank > 0) {
        if (!"vcov" %in% names(object)){
            constrain <- object$constrain
            eliminate <- object$eliminate
            needToElim <- seq(sum(!constrain[seq(eliminate)])[eliminate > 0])
            X <- model.matrix(object)[, !constrain, drop = FALSE]
            Info <- crossprod(X, as.vector(object$weights) * X)
            if (sum(constrain) > 0) {
                cov.unscaled <- array(0, dim = rep(length(constrain), 2),
                                      dimnames = rep(list(names(coefs)), 2))
                cov.unscaled[!constrain, !constrain] <-
                    MPinv(Info, eliminate = needToElim, onlyNonElim = FALSE)
            }
            else
                cov.unscaled <- MPinv(Info, eliminate = needToElim,
                                      onlyNonElim = FALSE)
            attr(cov.unscaled, "rank") <- NULL
        }
        else cov.unscaled  <- object$vcov
        cov.scaled <- dispersion * cov.unscaled
        estimable <- checkEstimable(object, diag(length(coefs)), ...)
        estimable[is.na(estimable)] <- FALSE
        sterr <- sqrt(diag(cov.scaled))
        is.na(sterr[!estimable]) <- TRUE
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
        cov.unscaled <- cov.scaled <- matrix(, 0, 0)
    }
    df.f <- nrow(coef.table)
    ans <- c(object[c("call", "terms", "family", "deviance", "aic",
                      "df.residual", "iter")],
             list(deviance.resid = residuals(object, type = "deviance"), 
                  coefficients = coef.table, dispersion = dispersion,
                  df = c(object$rank, object$df.residual, df.f),
                  eliminate = object$eliminate, 
                  cov.unscaled = cov.unscaled, cov.scaled = cov.scaled))
    if (correlation & object$rank > 0) {
        dd <- sqrt(diag(cov.unscaled))
        ans$correlation <- cov.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.gnm"
    return(ans)
}
