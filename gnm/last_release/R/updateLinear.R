updateLinear <- function(which, theta, y, mu, eta, offset, weights, family,
                         modelTools, X, eliminate) {
    dmu <- family$mu.eta(eta)
    vmu <- family$variance(mu)
    w <- weights * dmu * dmu / vmu
    theta[which] <- 0
    offsetVarPredictors <- modelTools$varPredictors(theta)
    offset <- offset + modelTools$predictor(offsetVarPredictors)
    z <- eta - offset + (y - mu)/dmu
    if (is.null(eliminate))
        naToZero(lm.wfit(X[,which, drop = FALSE], z, w)$coef)
    else
        suppressWarnings(glm.fit.e(X[,which, drop = FALSE], z,
                                   weights = w, intercept = FALSE,
                                   eliminate = eliminate, coefonly = TRUE))
}

