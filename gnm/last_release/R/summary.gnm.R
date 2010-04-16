summary.gnm <- function (object, dispersion = NULL, correlation = FALSE,
                         symbolic.cor = FALSE, with.eliminate = FALSE, ...)
{
    est.disp <- (!object$family$family %in% c("poisson", "binomial") &&
                 is.null(dispersion) && object$df.residual > 0)
    coefs <- parameters(object)
    if (with.eliminate) coefs <- c(attr(coef(object), "eliminated"), coefs)
    if (object$rank > 0) {
        cov.scaled <- vcov(object, dispersion = dispersion,
                           with.eliminate = with.eliminate)
        ## non-eliminated par only
        estimable <- checkEstimable(object, ...)
        estimable[is.na(estimable)] <- FALSE
        if (is.matrix(cov.scaled))
            sterr <- sqrt(diag(cov.scaled))
        else
            sterr <- diag(cov.scaled)
        is.na(sterr[!estimable]) <- TRUE
        if (with.eliminate){
            ## will need to do in gnmFit/checkestimable to get rank of whole model right in first place
            ## then can do object$rank - sum(estimable) to get rank of elim
            ## check estimability of eliminated coefficients
            browser()
            X <- cbind(1, model.matrix(object)[,!is.na(coef(object))])
            nelim <- nlevels(object$eliminate)
            rank <- object$rank - nelim
            estimable2 <- tapply(1:nrow(X), object$eliminate,
                                 function(ind) rankMatrix(X[ind,]) == rankMatrix(X[ind, -1]) + 1)
            ## OR
            X <- cbind(0, model.matrix(object)[,!is.na(coef(object))])
            nelim <- nlevels(object$eliminate)
            rank <- object$rank - nelim
            estimable2 <- tapply(1:nrow(X), object$eliminate,
                                 function(ind) {
                                     X[ind, 1] <- 1
                                     check <- rankMatrix(X) == rank + 1
                                     X[ind, 1] <- 0
                                     check
                                     })
            sterr <- c(ifelse(estimable2, sqrt(attr(cov.scaled, "varElim")), NA), sterr)
        }
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
    non.elim <- seq(object$coef) + nlevels(object$eliminate)*with.eliminate
    elim <- seq(length.out = nlevels(object$eliminate)*with.eliminate)
    ans <- c(object[c("call", "ofInterest", "family", "deviance", "aic",
                      "df.residual", "iter")],
             list(deviance.resid = residuals(object, type = "deviance"),
                  coefficients = coef.table[non.elim,],
                  eliminated = coef.table[elim,],
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
