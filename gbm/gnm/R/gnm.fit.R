gnm.fit <- function(modelTools, y, constrain, family = poisson(),
                    weights = rep.int(1, length(y)),
                    offset = rep.int(0, length(y)), nObs = length(y),
                    start = NULL,
                    control = gnm.control(...), x, vcov,
                    eliminate = numeric(0)) {
##    eliminate <- 1:101 in the backpain example
##    Need to sort out how start and constrain arguments interact with
##    eliminate: the coef and vcov components of the object relate only
##    to non-eliminated parameters.
    conv <- FALSE
    attempt <- 1
    repeat {
        if (any(is.na(start))) {
            theta <- modelTools$start()
            theta[!is.na(start)] <- start[!is.na(start)]
            theta[constrain] <- 0
            oneAtATime <- {!modelTools$classIndex %in%
                           c("Linear", "plugInStart") & is.na(start)}
            for (iter in seq(length = control$startit * any(oneAtATime))) {
                for (i in seq(theta)[oneAtATime]) {
                    if (constrain[i]) break
                    factorList <- modelTools$factorList(theta)
                    eta <- offset + modelTools$predictor(factorList)
                    mu <- family$linkinv(eta)
                    X <- modelTools$localDesignFunction(theta, factorList)
                    dmu <- family$mu.eta(eta)
                    vmu <- family$variance(mu)
                    w <- weights * dmu * dmu / vmu
                    Xi <- X[,i]
                    score <- crossprod((y - mu)/dmu, w * Xi)
                    gradient <- crossprod(w, Xi^2)
                    theta[i] <- as.vector(theta[i] + score/gradient)
                }
                if (control$trace){
                    dev <- sum(family$dev.resids(y, mu, weights))
                    cat("Startup iteration", iter,
                        ". Deviance = ", dev, "\n")
                }
            }
        }    
        else theta <- structure(ifelse(!constrain, start, 0),
                                names = names(modelTools$classIndex))
        for (iter in seq(control$maxit)) {
            factorList <- modelTools$factorList(theta)
            eta <- offset + modelTools$predictor(factorList)
            X <- modelTools$localDesignFunction(theta, factorList)
            mu <- family$linkinv(eta)
            dmu <- family$mu.eta(eta)
            vmu <- family$variance(mu)
            w <- weights * dmu * dmu / vmu
            if (control$trace){
                dev <- sum(family$dev.resids(y, mu, weights))
                cat("Iteration", iter, ". Deviance = ", dev, "\n")
            }
            if (is.nan(dev)) print("need to restart here!")
            z <- (y - mu)/dmu
            WX <- w * X
            score <- drop(crossprod(z, WX))
            diagInfo <- colSums(X * WX)
            if (all(abs(score) < control$epsilon*sqrt(diagInfo))){
                conv <- TRUE
                break
            }            
            Z <- cbind(z, X)
            WZ <- w * Z
            ZWZ <- crossprod(Z, WZ)
            ZWZinv <- try(gInvSymm(ZWZ, eliminate = 1 + eliminate,
                                   first.col.only = TRUE),
                          silent = TRUE)
            if (inherits(ZWZinv, "try-error")) break
            theChange <- - (ZWZinv[, 1] / ZWZinv[1, 1])[-1] 
            theta <- theta + theChange 
            theta[constrain] <- 0
        }
        if (!inherits(ZWZinv, "try-error") | all(!is.na(start))) break
        else {
            attempt <- attempt + 1
            if (attempt > 5) {
                cat("Fit attempted 5 times without success: terminating.\n")
                break
            }
            cat("Bad parameterisation: restarting.\n")
        }
    }
    if (!conv) {
        if (!all(is.finite(Info)))
            stop("Fit unsuccessful: values in information matrix are not ",
                 "all finite.\n")
        else
            warning("Fitting algorithm has either not converged or converged\n",
                    "to a non-solution of the likelihood equations.\n",
                    "Re-start gnm with coefficients of returned model.\n")
    }
    theta[constrain] <- NA
    WX <- w * X
    Info <- crossprod(X, WX)
    VCOV <- try(gInvSymm(Info, eliminate = eliminate, non.elim.only = TRUE),
                silent = TRUE)
    dev <- sum(family$dev.resids(y, mu, weights))
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs), mu, weights,
                                            dev) + 2 * attr(VCOV, "rank"))
    
    fit <- list(coefficients =
                if (length(eliminate) == 0) theta else theta[-eliminate],
                predictors = eta,
                fitted.values = mu,
                deviance = dev,
                aic = modelAIC,
                iter = iter,
                conv = conv,
                weights = w,
                residuals = z,
                df.residual = nObs - attr(VCOV, "rank"),
                rank = attr(VCOV, "rank"))
    if (x) fit$x <- structure(X, assign = modelTools$termAssign)
    if (vcov) {
        VCOV[constrain, constrain] <- 0
        fit$vcov <- VCOV
    }
    fit    
}

