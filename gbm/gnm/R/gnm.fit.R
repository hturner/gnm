gnm.fit <- function(modelTools, y, constrain, family = poisson(),
                     weights = rep.int(1, length(y)),
                     offset = rep.int(0, length(y)), nObs = length(y),
                     start = NULL, control = gnm.control(...), x, vcov) {
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
                    score <- crossprod(weights * (y - mu) * dmu / vmu, X[,i])
                    gradient <- crossprod(weights*dmu*dmu/vmu, (X[,i])^2)
                    theta[i] <- as.vector(theta[i] + score/gradient)
                }
                dev <- sum(family$dev.resids(y, mu, weights))
                if (control$trace)
                    cat("Startup iteration", iter,
                        ". Deviance = ", dev, "\n")    
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
            W <- diag(weights*dmu*dmu/vmu)
            WX <- crossprod(W, X)
            Info <- crossprod(X, WX)
            VCOV <- try(MPinv(Info, tol = 100*.Machine$double.eps),
                        silent = TRUE)
            if (inherits(VCOV, "try-error")) break
            dev <- sum(family$dev.resids(y, mu, weights))
            if (control$trace)
                cat("Iteration", iter, ". Deviance = ", dev, "\n")
            score <- drop(crossprod(weights*(y - mu)*dmu/vmu, X))
            if (all(abs(score) < control$epsilon*sqrt(diag(Info)))){
                conv <- TRUE
                break
            }
            theta <- theta + drop(crossprod(VCOV, score))
            theta[constrain] <- 0
        }
        if (!inherits(VCOV, "try-error") | all(!is.na(start))) break
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
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs), mu, weights,
                                            dev) + 2 * attr(VCOV, "rank"))
    fit <- list(coefficients = theta, predictors = eta, fitted.values = mu,
                deviance = dev, aic = modelAIC, iter = iter, conv = conv,
                weights = diag(W), residuals = (y - mu)/dmu,
                df.residual = nObs - attr(VCOV, "rank"),
                rank = attr(VCOV, "rank"))
    if (x) fit$x <- structure(X, assign = modelTools$termAssign)
    if (vcov) {
        VCOV[constrain, constrain] <- 0
        fit$vcov <- VCOV
    }
    fit    
}

