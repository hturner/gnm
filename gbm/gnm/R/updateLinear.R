updateLinear <- function(which, theta, y, mu, eta, offset, weights, family,
                         modelTools, X) {
    dmu <- family$mu.eta(eta)
    vmu <- family$variance(mu)
    w <- weights * dmu * dmu / vmu
    theta[which] <- NA
    offsetFactorList <- modelTools$factorList(theta, term = TRUE)
    offsetFactorList <- lapply(offsetFactorList, naToZero)
    offset <- offset + modelTools$predictor(offsetFactorList)
    z <- eta - offset + (y - mu)/dmu
    theta[which] <- suppressWarnings(naToZero(lm.wfit(X[,which, drop = FALSE],
                                                      z, w)$coef))
    theta
}
            
