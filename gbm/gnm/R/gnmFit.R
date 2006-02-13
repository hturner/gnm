"gnmFit" <-
    function (modelTools, y, constrain = FALSE, eliminate = 0,
              family = poisson(),
              weights = rep.int(1, length(y)),
              offset = rep.int(0, length(y)),
              nObs = length(y),
              start = rep.int(NA, length(y)),
              tolerance = 1e-4,
              iterStart = 2,
              iterMax = 500,
              trace = FALSE,
              verbose = FALSE,
              x = FALSE,
              termPredictors = FALSE,
              lsMethod = "svd")
{
    if (!(lsMethod %in% c("chol", "qr", "qr2", "svd"))) stop(
                "lsMethod must be one of chol, qr, svd")
    eps <- 100*.Machine$double.eps
    attempt <- 1
    dev <- numeric(2)
    if (verbose)
        width <- as.numeric(options("width"))
    repeat {
        status <- "not.converged"
        if (any(is.na(start))) {
            if (verbose == TRUE)
                prattle("Initialising", "\n", sep = "")
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
                theta[unspecifiedLin] <- quick.glm.fit(X[, unspecifiedLin], y,
                                                    weights = weights,
                                                    offset = offsetSpecified,
                                                    family = family,
                                                    eliminate = eliminate)
                constrain[is.na(theta)] <- TRUE
                theta <- naToZero(theta)
            }
            factorList <- modelTools$factorList(theta)
            eta <- offset + modelTools$predictor(factorList)
            mu <- family$linkinv(eta)
            dev[1] <- sum(family$dev.resids(y, mu, weights))
            if (trace)
                prattle("Initial Deviance = ", dev[1], "\n", sep = "")
            oneAtATime <- !linear & !specified
            for (iter in seq(length = iterStart * any(oneAtATime))) {
                if (verbose) {
                    if (iter == 1)
                        prattle("Running start-up iterations", "\n"[trace],
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
                    score <- crossprod(ifelse(abs(y - mu) < eps,
                                              0,
                                              (y - mu)/dmu),
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
                if (status == "not.converged" && any(linear)) {
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
                if (trace)
                    prattle("Start-up iteration ", iter, ". Deviance = ",
                            dev[1], "\n", sep = "")
                else if (verbose)
                    prattle(".")
                if (status == "bad.param")
                    break
                cat("\n"[iter == iterStart & verbose & !trace])
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
            if (trace)
                prattle("Initial Deviance = ", dev, "\n", sep = "")
        }
        if (status == "not.converged") {
            needToElim <- seq(sum(!constrain[seq(eliminate)])[eliminate > 0])
            X <-  modelTools$localDesignFunction(theta, factorList)
            rankX <- qr(X)$rank
            for (iter in seq(iterMax)) {
                if (verbose) {
                    if (iter == 1)
                        prattle("Running main iterations", "\n"[trace],
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
                wSqrt <- sqrt(w)
                W.X <- wSqrt * X
                score <- drop(crossprod(z, WX))
                diagInfo <- colSums(X * WX)
                if (diagInfo < 1e-20 ||
                    all(abs(score) < tolerance * sqrt(diagInfo))) {
                    status <- "converged"
                    break
                }
                if (iter > 1 && abs(diff(dev)) < 1e-16) {
                    status <- "stuck"
                    break
                }
                znorm <- sqrt(mean(z * z))#
                zscaled <- z/znorm#
                w.z <- wSqrt * zscaled#
                if (lsMethod %in% c("svd", "chol")) {
                    W.Z <- cbind(w.z, W.X)
                    ZWZ <- crossprod(W.Z)
                    ZWZinv <- MPinv((ZWZ),
                                    eliminate = 1 + needToElim,
                                    onlyFirstCol = TRUE,
                                    method = lsMethod)
                    theChange <- -(ZWZinv[, 1]/ZWZinv[1, 1])[-1] *
                  znorm#
                }
                browser()
                if (lsMethod == "qr") {
                    XWX <- crossprod(W.X)
                    theChange <- naToZero(qrSolve(XWX,
                                     crossprod(W.X, w.z))) * znorm#
                }
                if (lsMethod == "qr2") {
                    XWX <- crossprod(W.X)
                    theChange <- naToZero(qrSolve2(XWX,
                                     crossprod(W.X, w.z))) * znorm#
                }
                dev[2] <- dev[1]
                j <- 1
                while (dev[1] >= dev[2] && j < 11) {
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
                if (trace){
                    prattle("Iteration ", iter, ". Deviance = ", dev[1],
                            "\n", sep = "")
                }
                else if (verbose)
                    prattle(".")
                theta <- nextTheta
            }
        }
        if (status %in% c("converged", "not.converged")) {
            if (verbose)
                prattle("\n"[!trace], "Done\n", sep = "")
            break
        }
        else {
            if (verbose)
                message("\n"[!trace],
                        switch(status,
                               bad.param = "Bad parameterisation",
                               eta.not.finite =
                                 "Predictors are not all finite",
                               w.not.finite =
                                 "Iterative weights are not all finite",
                               no.deviance = "Deviance is NaN",
                               stuck = "Iterations are not converging"))
            attempt <- attempt + 1
            if (attempt > 5 || all(!is.na(start)) || modelTools$classID %in%
                                  c("Linear", "plugInStart"))
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
    Info <- crossprod(W.X)
    VCOV <- MPinv(Info, eliminate = needToElim, onlyNonElim = FALSE,
                  method = "svd")
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs),
                                            mu, weights, dev[1])
                                 + 2 * attr(VCOV, "rank"))
    fit <- list(coefficients = theta, constrain = constrain, residuals = z,
                fitted.values = mu, rank = attr(VCOV, "rank"), family = family,
                predictors = eta, deviance = dev[1], aic = modelAIC,
                iter = iter - 1, weights = w, prior.weights = weights,
                df.residual = nObs - sum(weights == 0) - attr(VCOV,"rank"),
                y = y, converged = status == "converged")
    if (x) {
        if (sum(constrain) > 0) {
            fit$x <- array(0, dim = c(nrow(X), length(theta)),
                           dimnames = list(NULL, names(theta)))
            fit$x[, !constrain] <- X
            attr(fit$x, "assign") <- modelTools$termAssign
        }
        else
            fit$x <- structure(X, assign = modelTools$termAssign)
    }
    if (termPredictors) {
        theta[is.na(theta)] <- 0
        factorList <- modelTools$factorList(theta, term = TRUE)
        fit$termPredictors <- modelTools$predictor(factorList, term = TRUE)
    }
    fit
}
