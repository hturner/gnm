"gnmFit" <-
    function (modelTools, y, constrain,
              family = poisson(),
              weights = rep.int(1, length(y)),
              offset = rep.int(0, length(y)),
              nObs = length(y),
              start,
              control = gnmControl(...),
              verbose = FALSE,
              x = FALSE,
              vcov = FALSE,
              termPredictors = FALSE)
{
    attempt <- 1
    if (verbose)
        width <- as.numeric(options("width"))
    repeat {
        status <- "not.converged"
        dev <- numeric(2)
        if (any(is.na(start))) {
            theta <- modelTools$start()
            theta[!is.na(start)] <- start[!is.na(start)]
            theta[constrain] <- 0
            modelTools$classID[is.na(theta)] <- "Linear"
            linear <- modelTools$classID == "Linear"
            specified <- !is.na(start) | modelTools$classID ==
                "plugInStart"
            unspecifiedLin <- linear & !specified & !constrain
            if (any(unspecifiedLin)) {
                thetaOffset <- theta
                thetaOffset[!specified] <- 0
                factorList <- modelTools$factorList(thetaOffset)
                offsetSpecified <- offset + modelTools$predictor(factorList)
                X <- modelTools$localDesignFunction(thetaOffset,
                                                    factorList)
                theta[unspecifiedLin] <- suppressWarnings(naToZero(
                                             glm.fit(X[, unspecifiedLin], y,
                                                     weights = weights,
                                                     offset = offsetSpecified,
                                                     family = family)$coef))
            }
            oneAtATime <- !linear & !specified & !constrain
            for (iter in seq(length = control$iterStart * any(oneAtATime))) {
                if (verbose) {
                    if (iter == 1)
                        cat("Running start-up iterations", "\n"[control$trace],
                            sep = "")
                    if ((iter + 25)%%width == (width - 1))
                        cat("\n")
                    if (!control$trace)
                        cat(".")
                }
                for (i in rep(seq(theta)[oneAtATime], 2)) {
                    factorList <- modelTools$factorList(theta)
                    eta <- offset + modelTools$predictor(factorList)
                    mu <- family$linkinv(eta)
                    X <- modelTools$localDesignFunction(theta,
                                                        factorList)
                    dmu <- family$mu.eta(eta)
                    vmu <- family$variance(mu)
                    w <- weights * dmu * dmu/vmu
                    Xi <- X[, i]
                    score <- crossprod((y - mu)/dmu, w * Xi)
                    gradient <- crossprod(w, Xi^2)
                    theta[i] <- as.vector(theta[i] + score/gradient)
                    if (!is.finite(theta[i])) {
                        status <- "bad.param"
                        break
                    }
                }
                if (status == "not.converged" & any(linear))
                    theta <- updateLinear(linear & !constrain,
                                          theta, y, offset, weights,
                                          family, modelTools,
                                          X)
                if (control$trace) {
                    dev <- sum(family$dev.resids(y, mu, weights))
                    cat("Start-up iteration ", iter, ". Deviance = ",
                        dev, "\n", sep = "")
                }
                if (status == "bad.param")
                    break
                cat("\n"[iter == control$iterStart & verbose &
                         !control$trace])
            }
        }
        else theta <- structure(ifelse(!constrain, start, 0),
                                names = names(modelTools$classID))
        for (iter in seq(control$iterMax)[status == "not.converged"]) {
            if (verbose) {
                if (iter == 1)
                    cat("Running main iterations", "\n"[control$trace],
                        sep = "")
                if ((iter + 21)%%width == (width - 1))
                    cat("\n")
                if (!control$trace)
                    cat(".")
            }
            factorList <- modelTools$factorList(theta)
            eta <- offset + modelTools$predictor(factorList)
            if (any(!is.finite(eta))) {
                status <- "eta.not.finite"
                break
            }
            X <- modelTools$localDesignFunction(theta, factorList)
            mu <- family$linkinv(eta)
            dmu <- family$mu.eta(eta)
            vmu <- family$variance(mu)
            w <- weights * dmu * dmu/vmu
            if (any(!is.finite(w))) {
                status <- "w.not.finite"
                break
            }
            dev[2] <- dev[1]
            dev[1] <- sum(family$dev.resids(y, mu, weights))
            if (control$trace)
                cat("Iteration ", iter, ". Deviance = ", dev[1],
                    "\n", sep = "")
            if (is.nan(dev[1])) {
                status <- "no.deviance"
                break
            }
            z <- (y - mu)/dmu
            WX <- w * X
            score <- drop(crossprod(z, WX))
            diagInfo <- colSums(X * WX)
            if (all(abs(score) < control$tolerance * sqrt(diagInfo) |
                    diagInfo < 1e-20)) {
                status <- "converged"
                break
            }
            if (iter > 1 & abs(diff(dev)) < 1e-16) {
                status <- "stuck"
                break
            }
            znorm <- sqrt(mean(z*z))
            zscaled <- z/znorm
            Z <- cbind(zscaled, X)
            WZ <- w * Z
            ZWZ <- crossprod(Z, WZ)
            ZWZinv <- MPinv(ZWZ,
                            eliminate = 1 + modelTools$eliminate,
                            onlyFirstCol = TRUE)
            theChange <- -(ZWZinv[, 1]/ZWZinv[1, 1])[-1] * znorm
            theta <- theta + theChange
            theta[constrain] <- 0
        }
        if (status %in% c("converged", "not.converged")) {
            if (verbose)
                cat("\n"[!control$trace], "Done\n", sep = "")
            break
        }
        else {
            if (verbose)
                message("\n"[!control$trace],
                        switch(status,
                               bad.param = "Bad parameterisation",
                               eta.not.finite = "Predictors are not all finite",
                               w.not.finite =
                               "Iterative weights are not all finite",
                               no.deviance = "Deviance is NaN",
                               stuck = "Iterations are not converging"))
            attempt <- attempt + 1
            if (attempt > 5 | all(!is.na(start) | modelTools$classID %in%
                                  c("Linear", "plugInStart")))
                return()
            else if (verbose)
                message("Restarting")
        }
    }
    if (status == "not.converged")
        warning("fitting algorithm has either not converged or converged\n",
                "to a non-solution of the likelihood equations: re-start \n",
                "gnm with coefficients of returned model\n")
    theta[constrain] <- NA
    if (exists("WX"))
        Info <- crossprod(X, WX)
    VCOV <- try(MPinv(Info, eliminate = modelTools$eliminate,
                      onlyNonElim = TRUE), silent = TRUE)
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs),
                                            mu, weights, dev[1])
                                 + 2 * attr(VCOV, "rank"))
    fit <- list(coefficients = theta, eliminate = length(modelTools$eliminate),
                predictors = eta, fitted.values = mu, deviance = dev[1],
                aic = modelAIC, iter = iter, conv = status == "converged",
                weights = w, residuals = z,
                df.residual = nObs - attr(VCOV,"rank"),
                rank = attr(VCOV, "rank"))
    if (x)
        fit$x <- structure(X, assign = modelTools$termAssign)
    if (vcov) {
        if (length(modelTools$eliminate))
            constrain <- constrain[-modelTools$eliminate]
        VCOV[constrain, constrain] <- 0
        fit$vcov <- VCOV
    }
    if (termPredictors) {
        theta[is.na(theta)] <- 0
        factorList <- modelTools$factorList(theta, term = TRUE)
        fit$termPredictors <- modelTools$predictor(factorList,
                                                   term = TRUE)
    }
    fit
}
