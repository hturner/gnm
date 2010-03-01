## returns vcov for the non-eliminated parameters
## use.eliminate now means whether or not to return extra components
## related to the eliminated parameters
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
    constrain <- object$constrain
    eliminate <- object$eliminate
    nelim <- nlevels(eliminate)
    w <- as.vector(object$weights)
    X <- model.matrix(object)
    ind <- !(seq_len(ncol(X)) %in% constrain)
    cov.unscaled <- array(0, dim = rep(ncol(X), 2),
                          dimnames = list(colnames(X), colnames(X)))
    if (!length(ind)) {
        if (nelim && use.eliminate) {
            Ti <- 1/sapply(split(w, eliminate), sum)
            attr(cov.unscaled, "varElim") <- dispersion * Ti
        }
        return(structure(cov.unscaled, dispersion = dispersion,
                         ofInterest = NULL, class = "vcov.gnm"))
    }

    if (length(constrain)) X <- X[, -constrain, drop = FALSE]
    W.X <- sqrt(w) * X
    if (object$rank == ncol(W.X)) {
        cov.unscaled[ind, ind] <- chol2inv(chol(crossprod(W.X)))
    } else {
        if (is.null(eliminate)) {
            cov.unscaled[ind, ind] <- MPinv(crossprod(W.X), method = "chol",
                                            rank = object$rank)
        } else {
            ## try without ridge and generalized inverse of Q
            Ti <- 1/sapply(split(w, eliminate), sum)
            U <- rowsum(sqrt(w) * W.X, eliminate)
            W <- crossprod(W.X)
            Ti.U <- Ti * U
            UTU <- crossprod(U, Ti.U)
            cov.unscaled[ind, ind] <- MPinv(W - UTU, method = "chol",
                                            rank = object$rank - nelim)
            if (use.eliminate) {
                rownames(Ti.U) <- names(object$elim.coefs)
                attr(cov.unscaled, "covElim") <- dispersion *
                    -Ti.U %*% cov.unscaled[ind, ind]
                attr(cov.unscaled, "varElim") <- dispersion *
                    -diag(tcrossprod(attr(cov.unscaled, "covElim"), Ti.U)) + Ti
            }
        }
    }
    structure(dispersion * cov.unscaled, dispersion = dispersion,
              ofInterest = ofInterest(object), class = "vcov.gnm")
}
