"gnmFit" <-
    function (modelTools, y,
              constrain = numeric(0),
              constrainTo = numeric(length(constrain)),
              eliminate = 0,
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
              lsMethod = "qr",
              ridge = 1e-8)
{
    if (!(lsMethod %in% c("chol", "qr"))) stop(
                "lsMethod must be chol or qr")
    eps <- 100*.Machine$double.eps
    attempt <- 1
    dev <- numeric(2)
    if (verbose)
        width <- as.numeric(options("width"))
    isConstrained <- is.element(seq(start), constrain)
    XWX <- NULL
    repeat {
        status <- "not.converged"
        if (any(is.na(start))) {
            if (verbose == TRUE)
                prattle("Initialising", "\n", sep = "")
            theta <- modelTools$start()
            theta[!is.na(start)] <- start[!is.na(start)]
            theta[constrain] <- constrainTo
            modelTools$classID[is.na(theta)] <- "Linear"
            linear <- modelTools$classID == "Linear"
            specified <- !is.na(start) | modelTools$classID ==
                "plugInStart" | isConstrained
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
                if (sum(is.na(theta)) > length(constrain)) {
                    extra <- setdiff(which(is.na(theta)), constrain)
                    isConstrained[extra] <- TRUE
                    ind <- order(c(constrain, extra))
                    constrain <- c(constrain, extra)[ind]
                    constrainTo <- c(constrainTo, numeric(length(extra)))[ind]
                }
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
                        which <- seq(theta)[linear & !isConstrained]
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
            theta <- structure(replace(start, constrain, constrainTo),
                               names = names(modelTools$classID))
            factorList <- modelTools$factorList(theta)
            eta <- offset + modelTools$predictor(factorList)
            if (any(!is.finite(eta))) {
                stop("Values of 'start' and 'constrain' produce non-finite ",
                     "predictor values")
            }
            mu <- family$linkinv(eta)
            dev[1] <- sum(family$dev.resids(y, mu, weights))
            if (trace)
                prattle("Initial Deviance = ", dev, "\n", sep = "")
        }
        if (status == "not.converged") {
            needToElim <- seq(length.out = eliminate)
            X <-  modelTools$localDesignFunction(theta, factorList)
            X <- X[, !isConstrained, drop = FALSE]
            pns <- rep(nrow(X), ncol(X))
            ridge <- c(0, rep(ridge, ncol(X)))
            for (iter in seq(iterMax)) {
                if (any(!is.finite(X))){
                    status <- "X.not.finite"
                    break
                }
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

                wSqrt <- sqrt(w)
                W.X <- wSqrt * X
                w.z <- wSqrt * z
                score <- drop(crossprod(w.z, W.X))
                diagInfo <- colSums(W.X * W.X)
                Xscales <- pmax(1e-3, sqrt(diagInfo)) ## to allow for zeros
                W.X.scaled <- W.X / rep(Xscales, pns)
                if (all(diagInfo < 1e-20) ||
                    all(abs(score) <
                        tolerance * sqrt(tolerance + diagInfo))) {
                    status <- "converged"
                    break
                }
                znorm <- sqrt(sum(w.z * w.z))
                w.z <- w.z/znorm
                if (lsMethod == "chol") {
                    W.Z <- cbind(w.z, W.X.scaled)
                    if (eliminate > 0){
                        Tvec <- rep(1, eliminate) + ridge[1 + needToElim]
                        Wmat <- W.Z[, -(1 + needToElim), drop = FALSE]
                        Umat <- crossprod(W.Z[, 1 + needToElim, drop = FALSE],
                                          Wmat)
                        Wmat <- crossprod(Wmat)
                        diag(Wmat) <- diag(Wmat) + ridge[-(1 + needToElim)]
                        ZWZinv <- cholInv1(Wmat, Tvec, Umat, 1 + needToElim)
                    } else {
                        ZWZ <- crossprod(W.Z)
                        theDiagonal <- diag(ZWZ)
                        diag(ZWZ) <- theDiagonal + ridge
                        ZWZinv <- cholInv1(ZWZ)
                    }
                    theChange <- -ZWZinv[-1]/ZWZinv[1] * znorm / Xscales
                } else { ## lsMethod is "qr"
                    XWX <- crossprod(W.X.scaled)
                    theDiagonal <- diag(XWX)
                    XWXridge <- XWX
                    diag(XWXridge) <- theDiagonal + ridge[-1]
                    Qr <- qr(XWXridge, tol = 1e-20)
                    theChange <- solve(Qr, crossprod(W.X.scaled, w.z)) *
                        znorm / Xscales
                }
                dev[2] <- dev[1]
                j <- 1
                while (dev[1] >= dev[2] && j < 11) {
                    nextTheta <- replace(theta, !isConstrained,
                                         theta[!isConstrained] + theChange)
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
                X <- modelTools$localDesignFunction(theta, factorList)
                X <- X[, !isConstrained, drop = FALSE]
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
                               X.not.finite =
                                 "Local design matrix has infinite elements",
                               no.deviance = "Deviance is NaN"))
            attempt <- attempt + 1
            if (attempt > 5 || all(!is.na(start)) || modelTools$classID %in%
                                  c("Linear", "plugInStart"))
                return()
            else if (verbose)
                message("Restarting")
        }
    }
    theta[constrain] <- NA
    sv.tolerance <-  100 * .Machine$double.eps
    if (is.null(XWX)) XWX <- crossprod(W.X.scaled)
    Svd <- svd(XWX, nu = 0, nv = 0)
    theRank <- sum(Svd$d > max(sv.tolerance * Svd$d[1], 0))
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs),
                                            mu, weights, dev[1])
                                 + 2 * theRank)
    fit <- list(coefficients = theta, constrain = constrain,
                constrainTo = constrainTo, residuals = z, fitted.values = mu,
                rank = theRank, family = family, predictors = eta,
                deviance = dev[1], aic = modelAIC,
                iter = iter - (iter != iterMax), weights = w,
                prior.weights = weights,
                df.residual = nObs - theRank,
                y = y)
    if (status == "not.converged") {
        warning("Fitting algorithm has either not converged or converged\n",
                "to a non-solution of the likelihood equations.\n",
                "Use exitInfo() for numerical details of last iteration.\n")
        fit$converged <- structure(FALSE, score = score, criterion =
                                   tolerance * sqrt(tolerance + diagInfo))
    }
    else
        fit$converged <- TRUE

    if (x) {
        if (length(constrain) > 0) {
            fit$x <- array(0, dim = c(nrow(X), length(theta)),
                           dimnames = list(NULL, names(theta)))
            fit$x[, !isConstrained] <- X
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
