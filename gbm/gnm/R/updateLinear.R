updateLinear <- function(which, theta, y, offset, weights, family, modelTools,
                         X = NULL) {
    factorList <- modelTools$factorList(theta)
    eta <- offset + modelTools$predictor(factorList)
    mu <- family$linkinv(eta)
    if (is.null(X))
        X <- modelTools$localDesignFunction(theta, factorList)
    dmu <- family$mu.eta(eta)
    vmu <- family$variance(mu)
    w <- weights * dmu * dmu / vmu
    theta[which] <- 0
    remainderFactorList <- modelTools$factorList(theta)
    offsetRemainder <- offset +
        modelTools$predictor(remainderFactorList)
    z <- eta - offsetRemainder + (y - mu)/dmu
    theta[which] <- suppressWarnings(naToZero(lm.wfit(X[,which], z, w)$coef))
    theta
}
            
