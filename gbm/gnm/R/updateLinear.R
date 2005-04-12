updateLinear <- function(which, theta, y, offset, weights, family,
                         modelTools, X) {
    factorList <- modelTools$factorList(theta)
    eta <- offset + modelTools$predictor(factorList)
    mu <- family$linkinv(eta)
    dmu <- family$mu.eta(eta)
    vmu <- family$variance(mu)
    w <- weights * dmu * dmu / vmu
    theta[which] <- 0
    offsetFactorList <- modelTools$factorList(theta)
    offset <- offset + modelTools$predictor(offsetFactorList)
    z <- eta - offset + (y - mu)/dmu
    theta[which] <- suppressWarnings(naToZero(lm.wfit(X[,which], z, w)$coef))
    theta
}
            
