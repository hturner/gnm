gnm.fit <- function(modelTools, y, constrain, family = poisson(),
                    weights = rep.int(1, length(y)),
                    offset = rep.int(0, length(y)), nObs = length(y),
                    start = NULL,
                    control = gnm.control(...), x, vcov, term.predictors) {
##    eliminate <- 1:101 in the backpain example
##    Need to sort out how start and constrain arguments interact with
##    eliminate: the coef and vcov components of the object relate only
##    to non-eliminated parameters.
    attempt <- 1
    repeat {
        status <- "not.converged"
        dev <- numeric(2)
        if (any(is.na(start))) {
            theta <- modelTools$start()
            theta[!is.na(start)] <- start[!is.na(start)]
            theta[constrain] <- 0
            linear <- modelTools$classID == "Linear"
            freeLin <- linear & is.na(start)
            theta <- updateLinear(freeLin, theta, y, offset, weights, family,
                                  modelTools)
            oneAtATime <- {!linear & modelTools$classID != "plugInStart" &
                          is.na(start)}
            for (iter in seq(length = control$startit * any(oneAtATime))) {
                for (i in rep(seq(theta)[oneAtATime], 2)) {
                    if (constrain[i]) break
                    factorList <- modelTools$factorList(theta)
                    eta <- offset + modelTools$predictor(factorList)
                    mu <- family$linkinv(eta)
                    X <- modelTools$localDesignFunction(theta, factorList)
                    dmu <- family$mu.eta(eta)
                    vmu <- family$variance(mu)
                    w <- weights * dmu * dmu / vmu
                    Xi <- X[, i]
                    score <- crossprod((y - mu)/dmu, w * Xi)
                    gradient <- crossprod(w, Xi^2)
                    theta[i] <- as.vector(theta[i] + score/gradient)
                    if (!is.finite(theta[i])) {
                        status <- "bad.param"
                        break
                    } 
                }
                if (status == "not.converged") 
                    theta <- updateLinear(linear, theta, y, offset, weights,
                                          family, modelTools, X)
                if (control$trace){
                    dev <- sum(family$dev.resids(y, mu, weights))
                    cat("Startup iteration", iter,
                        ". Deviance = ", dev, "\n")
                }
                if (status == "bad.param") break
            }
        }    
        else theta <- structure(ifelse(!constrain, start, 0),
                                names = names(modelTools$classID))
        for (iter in seq(control$maxit)[status == "not.converged"]) {
            factorList <- modelTools$factorList(theta)
            eta <- offset + modelTools$predictor(factorList)
            X <- modelTools$localDesignFunction(theta, factorList)
            mu <- family$linkinv(eta)
            dmu <- family$mu.eta(eta)
            vmu <- family$variance(mu)
            w <- weights * dmu * dmu / vmu
            if (any(!is.finite(w))) {
                status <- "not.finite"
                break
            }
            dev[2] <- dev[1]
            dev[1] <- sum(family$dev.resids(y, mu, weights))
            if (control$trace)
                cat("Iteration", iter, ". Deviance = ", dev[1], "\n")
            if (is.nan(dev[1])) {
                status <- "no.deviance"
                break
            }
            z <- (y - mu)/dmu
            WX <- w * X
            score <- drop(crossprod(z, WX))
            diagInfo <- colSums(X * WX)
            if (all(abs(score) < control$epsilon*sqrt(diagInfo) |
                    diagInfo < 1e-20)){
                status <- "converged"
                break
            }
            if (iter >1 & abs(diff(dev)) < 1e-16) {
                status <- "stuck"
                break
            }
            Z <- cbind(z, X)
            WZ <- w * Z
            ZWZ <- crossprod(Z, WZ)
            ZWZinv <- gInvSymm(ZWZ, eliminate = 1 + modelTools$eliminate,
                                   first.col.only = TRUE)
            theChange <- - (ZWZinv[, 1] / ZWZinv[1, 1])[-1] 
            theta <- theta + theChange 
            theta[constrain] <- 0
        }
        if (status %in% c("converged", "not.converged") | all(!is.na(start)))
            break
        else {
            attempt <- attempt + 1
            cat(switch(status,
                       "bad.param" = "Bad parameterisation",
                       "not.finite" = "Iterative weights are not all finite",
                       "no.deviance" = "Deviance is NaN",
                       "stuck" = "Iterations are not converging"))
            if (attempt > 5) {
                cat(".\nFit attempted 5 times without success: terminating.\n")
                break
            }
            else cat(": restarting. \n")
            #cat("Bad parameterisation: restarting.\n")
        }
    }
    if (status == "not.converged")
        warning("Fitting algorithm has either not converged or converged\n",
                "to a non-solution of the likelihood equations.\n",
                "Re-start gnm with coefficients of returned model.\n")
    theta[constrain] <- NA
    if (exists("WX")) Info <- crossprod(X, WX)
    VCOV <- try(gInvSymm(Info, eliminate = modelTools$eliminate,
                         non.elim.only = TRUE), silent = TRUE)
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs), mu, weights,
                                            dev[1]) + 2 * attr(VCOV, "rank"))
    
    fit <- list(coefficients = if (length(modelTools$eliminate) == 0) theta
                else theta[-modelTools$eliminate],
                predictors = eta,
                fitted.values = mu,
                deviance = dev[1],
                aic = modelAIC,
                iter = iter,
                conv = status == "converged",
                weights = w,
                residuals = z,
                df.residual = nObs - attr(VCOV, "rank"),
                rank = attr(VCOV, "rank"))
    if (x) fit$x <- structure(X, assign = modelTools$termAssign)
    if (vcov) {
        VCOV[constrain, constrain] <- 0
        fit$vcov <- VCOV
    }
    if (term.predictors) {
        factorList <- modelTools$factorList(theta, term = TRUE)
        fit$term.predictors <- modelTools$predictor(factorList, term = TRUE)
    }
    fit    
}

