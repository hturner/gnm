gnm.fit <- function(modelTools, y, constrain, family = poisson(),
                    weights = rep.int(1, length(y)),
                    offset = rep.int(0, length(y)), nObs = length(y),
                    start = NULL,
                    control = gnm.control(...), x, vcov, term.predictors,
                    eliminate = numeric(0)) {
##    eliminate <- 1:101 in the backpain example
##    Need to sort out how start and constrain arguments interact with
##    eliminate: the coef and vcov components of the object relate only
##    to non-eliminated parameters.
    attempt <- 1
    eliminateChecked <- length(eliminate) == 0
    repeat {
        status <- "not.converged"
        dev <- numeric(2)
        if (any(is.na(start))) {
            theta <- modelTools$start()
            theta[!is.na(start)] <- start[!is.na(start)]
            theta[constrain] <- 0
            linear <- modelTools$classIndex == "Linear"
            # update linear terms, offsetting any Nonlin/fixed
            #factorList <- modelTools$factorList(theta)
#            eta <- offset + modelTools$predictor(factorList)
#            mu <- family$linkinv(eta)
#            X <- modelTools$localDesignFunction(theta, factorList)
#            if (!eliminateChecked) {
#                Xelim <- crossprod(X[, eliminate])
#                if (any(abs(Xelim[lower.tri(Xelim)])) > 1e-15) stop(
#                "eliminated parameters must correspond to levels of a factor")
#                eliminateChecked <- TRUE
#            }
#            dmu <- family$mu.eta(eta)
#            vmu <- family$variance(mu)
#            w <- weights * dmu * dmu / vmu
#            freeLin <- linear & is.na(start)
#            theta[freeLin] <- 0
#            remainderFactorList <- modelTools$factorList(theta)
#            offsetRemainder <- offset +
#                modelTools$predictor(remainderFactorList)
#            z <- eta - offsetRemainder + (y - mu)/dmu
#            theta[freeLin] <-
#                suppressWarnings(naToZero(lm.wfit(X[,freeLin], z, w)$coef))
            # loop updating Nonlin and lin terms except ones provided
            oneAtATime <- {!linear & modelTools$classIndex != "plugInStart" &
                          is.na (start)}
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
                if (status == "not.converged") {
                    factorList <- modelTools$factorList(theta)
                    eta <- offset + modelTools$predictor(factorList)
                    mu <- family$linkinv(eta)
                    dmu <- family$mu.eta(eta)
                    vmu <- family$variance(mu)
                    w <- weights * dmu * dmu / vmu
                    theta[linear] <- 0
                    nonlinearFactorList <- modelTools$factorList(theta)
                    offsetNonlin <- offset +
                        modelTools$predictor(nonlinearFactorList)
                    z <- eta - offsetNonlin + (y - mu)/dmu
                    theta[linear] <-
                        suppressWarnings(naToZero(lm.wfit(X[,linear],
                                                          z, w)$coef))
                }
                if (control$trace){
                    dev <- sum(family$dev.resids(y, mu, weights))
                    cat("Startup iteration", iter,
                        ". Deviance = ", dev, "\n")
                }
                if (status == "bad.param") break
            }
        }    
        else theta <- structure(ifelse(!constrain, start, 0),
                                names = names(modelTools$classIndex))
        for (iter in seq(control$maxit)[status == "not.converged"]) {
            factorList <- modelTools$factorList(theta)
            eta <- offset + modelTools$predictor(factorList)
            X <- modelTools$localDesignFunction(theta, factorList)
            if (!eliminateChecked) {
                Xelim <- crossprod(X[, eliminate])
                if (any(abs(Xelim[lower.tri(Xelim)])) > 1e-15) stop(
                "eliminated parameters must correspond to levels of a factor")
                eliminateChecked <- TRUE
            }
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
            ZWZinv <- gInvSymm(ZWZ, eliminate = 1 + eliminate,
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
    VCOV <- try(gInvSymm(Info, eliminate = eliminate, non.elim.only = TRUE),
                silent = TRUE)
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs), mu, weights,
                                            dev[1]) + 2 * attr(VCOV, "rank"))
    
    fit <- list(coefficients =
                if (length(eliminate) == 0) theta else theta[-eliminate],
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

