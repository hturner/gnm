updateLinear <- function(which, theta, y, mu, eta, offset, weights, family,
                         modelTools, X) {
    dmu <- family$mu.eta(eta)
    vmu <- family$variance(mu)
    w <- weights * dmu * dmu / vmu
    theta[which] <- NA
    offsetVarPredictors <- modelTools$varPredictors(theta)
    offsetVarPredictors <- lapply(offsetVarPredictors, naToZero)
    offset <- offset + modelTools$predictor(offsetVarPredictors)
    z <- eta - offset + (y - mu)/dmu
    theta[which] <- suppressWarnings(naToZero(lm.wfit(X[,which, drop = FALSE],
                                                      z, w)$coef))
    theta
}
            
