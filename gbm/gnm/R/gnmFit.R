"gnmFit" <-
    function (modelTools, y, constrain, eliminate,
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
    eps <- 100*.Machine$double.eps
    attempt <- 1
    dev <- numeric(2)
    if (verbose)
        width <- as.numeric(options("width"))
    theta <- seq(start)
    X <- modelTools$localDesignFunction(theta, modelTools$factorList(theta))
    if (eliminate)
        if (nrow(unique(X[, seq(eliminate)])) > eliminate)
            stop("'eliminate' formula is not equivalent to single factor")
    unestimable <- apply(X == 0, 2, all)
    constrain[unestimable] <- TRUE
    repeat {
        status <- "not.converged"
        if (any(is.na(start))) {
            if (verbose == TRUE)
                cat("Initialising", "\n", sep = "")
            theta <- modelTools$start()
            theta[!is.na(start)] <- start[!is.na(start)]
            theta[constrain] <- 0
            modelTools$classID[is.na(theta)] <- "Linear"
            linear <- modelTools$classID == "Linear"
            specified <- !is.na(start) | modelTools$classID ==
                "plugInStart" | constrain
            unspecifiedLin <- seq(theta)[linear & !specified]
            if (any(unspecifiedLin)) {
                thetaOffset <- theta
                thetaOffset[!specified] <- NA
                factorList <- modelTools$factorList(thetaOffset,
                                                    term = TRUE)
                factorList <- lapply(factorList, naToZero)
                offsetSpecified <- offset + modelTools$predictor(factorList)
                X <- modelTools$localDesignFunction(thetaOffset,
                                                    factorList)
                theta[unspecifiedLin] <-
                    suppressWarnings(naToZero(glm.fit(X[, unspecifiedLin], y,
                                                      weights = weights,
                                                      offset = offsetSpecified,
                                                      family = family)$coef))
            }
            factorList <- modelTools$factorList(theta)
            eta <- offset + modelTools$predictor(factorList)
            mu <- family$linkinv(eta)
            dev[1] <- sum(family$dev.resids(y, mu, weights))
            if (control$trace)
                cat("Initial Deviance = ", dev[1], "\n", sep = "")
            oneAtATime <- !linear & !specified
            for (iter in seq(length = control$iterStart * any(oneAtATime))) {
                if (verbose) {
                    if (iter == 1)
                        cat("Running start-up iterations", "\n"[control$trace],
                            sep = "")
                    if ((iter + 25)%%width == (width - 1))
                        cat("\n")
                }
                for (i in rep(seq(theta)[oneAtATime], 2)) {
                    dmu <- family$mu.eta(eta)
                    vmu <- family$variance(mu)
                    w <- weights * ifelse(abs(dmu) < eps, 0, dmu * dmu/vmu)
                    Xi <- modelTools$localDesignFunction(theta,
                                                         factorList, i)
                    score <- crossprod(ifelse(abs(y - mu) < eps, 0, (y - mu)/dmu),
                                       w * Xi)
                    gradient <- crossprod(w, Xi^2)
                    theta[i] <- as.vector(theta[i] + score/gradient)
                    if (!is.finite(theta[i])) {
                        status <- "bad.param"
                        break
                    }
                    factorList <- modelTools$factorList(theta)
                    eta <- offset + modelTools$predictor(factorList)
                    mu <- family$linkinv(eta)
                }
                if (status == "not.converged" & any(linear)) {
                    if (iter == 1) {
                        which <- seq(theta)[linear & !constrain]
                        if(!exists("X"))
                            X <- modelTools$localDesignFunction(theta,
                                                                factorList)
                    }
                    theta <- updateLinear(which, theta, y, mu, eta, offset,
                                          weights, family, modelTools, X)
                    factorList <- modelTools$factorList(theta)
                    eta <- offset + modelTools$predictor(factorList)
                    mu <- family$linkinv(eta)
                }
                dev[1] <- sum(family$dev.resids(y, mu, weights))
                if (control$trace)
                    cat("Start-up iteration ", iter, ". Deviance = ",
                        dev[1], "\n", sep = "")
                else if (verbose)
                    cat(".")
                if (status == "bad.param")
                    break
                cat("\n"[iter == control$iterStart & verbose &
                         !control$trace])
            }
        }
        else {
            theta <- structure(ifelse(!constrain, start, 0),
                                names = names(modelTools$classID))
            factorList <- modelTools$factorList(theta)
            eta <- offset + modelTools$predictor(factorList)
            if (any(!is.finite(eta))) {
                status <- "eta.not.finite"
                break
            }
            mu <- family$linkinv(eta)
            dev[1] <- sum(family$dev.resids(y, mu, weights))
            if (control$trace)
                cat("Initial Deviance = ", dev, "\n", sep = "")
        }
        if (status == "not.converged") {
            needToElim <- seq(sum(!constrain[seq(eliminate)]))[eliminate > 0]
            for (iter in seq(control$iterMax)) {
                if (verbose) {
                    if (iter == 1)
                        cat("Running main iterations", "\n"[control$trace],
                            sep = "")
                    if ((iter + 21)%%width == (width - 1))
                        cat("\n")
                }
                dmu <- family$mu.eta(eta)
                z <- ifelse(abs(dmu) < eps, 0, (y - mu)/dmu)
                vmu <- family$variance(mu)
                w <- weights * ifelse(abs(dmu) < eps, 0, dmu * dmu/vmu)
                if (any(!is.finite(w))) {
                    status <- "w.not.finite"
                    break
                }
                X <- modelTools$localDesignFunction(theta, factorList)
                X <- X[, !constrain, drop = FALSE]
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
                                eliminate = 1 + needToElim,
                                onlyFirstCol = TRUE)
                theChange <- -(ZWZinv[, 1]/ZWZinv[1, 1])[-1] * znorm
                dev[2] <- dev[1]
                j <- 1
                while (dev[1] >= dev[2] & j < 11) {
                    nextTheta <- replace(theta, !constrain,
                                         theta[!constrain] + theChange)
                    factorList <- modelTools$factorList(nextTheta)
                    eta <- offset + modelTools$predictor(factorList)
                    if (any(!is.finite(eta))) {
                        status <- "eta.not.finite"
                        break
                    }
                    mu <- family$linkinv(eta)
                    dev[1] <- sum(family$dev.resids(y, mu, weights))
                    if (is.nan(dev[1])) {
                        status <- "no.deviance"
                        break
                    }
                    theChange <- theChange/2
                    j <- j + 1
                }
                if (control$trace){
                    cat("Iteration ", iter, ". Deviance = ", dev[1],
                        "\n", sep = "")
                }
                else if (verbose)
                    cat(".")
                theta <- nextTheta
            }
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
    Info <- crossprod(X, WX)
    VCOV <- MPinv(Info, eliminate = needToElim, onlyNonElim = TRUE)
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs),
                                            mu, weights, dev[1])
                                 + 2 * attr(VCOV, "rank"))
    fit <- list(coefficients = theta, eliminate = eliminate,
                constrain = constrain, predictors = eta, fitted.values = mu,
                deviance = dev[1], aic = modelAIC, iter = iter - 1,
                conv = status == "converged", weights = w, residuals = z,
                df.residual = nObs - sum(weights == 0) - attr(VCOV,"rank"),
                rank = attr(VCOV, "rank"))
    if (x) {
        if (sum(constrain) > 0) {
            fit$x <- array(0, dim = c(nrow(X), length(theta)),
                           dimnames = list(NULL, names(theta)))
            attr(fit$x, "assign") <- modelTools$termAssign
            fit$x[, !constrain] <- X
        }
        else
            fit$x <- structure(X, assign = modelTools$termAssign)
    }
    if (vcov) {
        if (eliminate)
            constrain <- constrain[-seq(eliminate)]
        if (sum(constrain) > 0) {
            fit$vcov <- array(0, dim = rep(length(theta), 2),
                              dimnames = rep(list(names(theta)), 2))
            fit$vcov[!constrain, !constrain] <- VCOV
        }
        else fit$vcov <- VCOV[, , drop = FALSE]
    }
    if (termPredictors) {
        theta[is.na(theta)] <- 0
        factorList <- modelTools$factorList(theta, term = TRUE)
        fit$termPredictors <- modelTools$predictor(factorList, term = TRUE)
    }
    fit
}
