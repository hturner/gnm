vcov.gnm <-  function(object, dispersion = NULL, use.eliminate = TRUE, ...){
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
    coefNames <- names(coef(object))
    constrain <- object$constrain
    eliminate <- object$eliminate
    needToElim <- seq(length.out = eliminate)
    isConstrained <- is.element(seq(coefNames), constrain)
    X <- model.matrix(object)[, !isConstrained, drop = FALSE]
    w <- as.vector(object$weights)
    W.X <- sqrt(w) * X
    cov.unscaled <- array(0, dim = rep(length(coefNames), 2),
                          dimnames = rep(list(coefNames), 2))
    if (object$rank == ncol(W.X)) {
        invInfo <- chol2inv(chol(crossprod(W.X)))
    } else {
        if (eliminate == 0 || !use.eliminate) {
            invInfo <- MPinv(crossprod(W.X), method = "chol",
                             rank = object$rank)
        } else {
            Tvec <- colSums(w * X[, seq(eliminate), drop = FALSE])
            Wmat <- W.X[, -needToElim, drop = FALSE]
            Umat <- crossprod(W.X[, needToElim, drop = FALSE], Wmat)
            Wmat <- crossprod(Wmat)
            invInfo <- MPinv(list(Wmat, Tvec, Umat), eliminate = needToElim,
                             method = "chol", rank = object$rank)
        }
    }
    cov.unscaled[!isConstrained, !isConstrained] <- invInfo
    attr(cov.unscaled, "rank") <- NULL
    structure(dispersion * cov.unscaled, dispersion = dispersion,
              ofInterest = ofInterest(object), class = "vcov.gnm")
}
