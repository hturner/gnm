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
    needToElim <- seq(length.out = eliminate - sum(constrain[,1] < eliminate))
    isConstrained <- is.element(seq(coef(object)), constrain[,1])
    X <- model.matrix(object)[, !isConstrained, drop = FALSE]
    Info <- crossprod(X, c(object$weights) * X)
    if (sum(constrain) > 0) {
        cov.unscaled <- array(0, dim = rep(length(isConstrained), 2),
                              dimnames = rep(list(names(coef(object))), 2))
        cov.unscaled[!isConstrained, !isConstrained] <-
            MPinv(Info, eliminate = needToElim, onlyNonElim = FALSE)
    }
    else
        cov.unscaled <- MPinv(Info, eliminate = needToElim, onlyNonElim = FALSE)
    attr(cov.unscaled, "rank") <- NULL
    structure(dispersion * cov.unscaled, dispersion = dispersion,
              eliminate = eliminate, class = "vcov.gnm")
}
