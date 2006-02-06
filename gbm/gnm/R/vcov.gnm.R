vcov.gnm <-  function(object, dispersion = NULL, ...){
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
    constrain <- object$constrain
    eliminate <- object$eliminate
    needToElim <- seq(sum(!constrain[seq(eliminate)])[eliminate > 0])
    X <- model.matrix(object)[, !constrain, drop = FALSE]
    Info <- crossprod(X, c(object$weights) * X)
    if (sum(constrain) > 0) {
        cov.unscaled <- array(0, dim = rep(length(constrain), 2),
                              dimnames = rep(list(names(coef(object))), 2))
        cov.unscaled[!constrain, !constrain] <-
            MPinv(Info, eliminate = needToElim, onlyNonElim = FALSE)
    }
    else
        cov.unscaled <- MPinv(Info, eliminate = needToElim, onlyNonElim = FALSE)
    attr(cov.unscaled, "rank") <- NULL
    structure(dispersion * cov.unscaled, dispersion = dispersion,
              eliminate = eliminate, class = "vcov.gnm")
}
