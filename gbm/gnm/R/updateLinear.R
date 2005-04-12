updateLinear <- function(which, add.to.offset, theta, y, offset, weights,
                         family, modelTools, X = NULL) {
    factorList <- modelTools$factorList(theta)
    eta <- offset + modelTools$predictor(factorList)
    mu <- family$linkinv(eta)
    if (is.null(X))
        X <- modelTools$localDesignFunction(theta, factorList)
    dmu <- family$mu.eta(eta)
    vmu <- family$variance(mu)
    w <- weights * dmu * dmu / vmu
    if (any(add.to.offset)) {
        theta.offset <- theta
        theta.offset[!add.to.offset] <- 0
        offsetFactorList <- modelTools$factorList(theta.offset)
        offset <- offset + modelTools$predictor(offsetFactorList)
    }
    z <- eta - offset + (y - mu)/dmu
    theta[which] <- suppressWarnings(naToZero(lm.wfit(X[,which], z, w)$coef))
    theta
}
            
