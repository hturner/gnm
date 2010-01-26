## simplified gnmFit
## ignore constrain, skip starting iterations (use start), assume nonlinear & method = chol
"gnmFit.e" <-
    function (modelTools, y,
              constrain = numeric(0),
              constrainTo = numeric(length(constrain)),
              eliminate = NULL, # now a factor
              family = poisson(),
              weights = rep.int(1, length(y)),
              offset = rep.int(0, length(y)),
              nObs = length(y),
              start = rep.int(NA, length(y)),
              etastart = NULL,
              mustart = NULL,
              tolerance = 1e-6,
              iterStart = 2,
              iterMax = 500,
              trace = FALSE,
              verbose = FALSE,
              x = FALSE,
              termPredictors = FALSE,
              ridge = 1e-8)
{
    eps <- 100*.Machine$double.eps
    attempt <- 1
    dev <- numeric(2)
    if (verbose)
        width <- as.numeric(options("width"))
    isConstrained <- is.element(seq(length(start)), constrain)
    XWX <- NULL
    repeat {
        status <- "not.converged"
        ## then need starts for elimCoef - use quick.glm.fit
        alpha <- rep.int(1, nlevels(eliminate))
        theta <- structure(replace(start, constrain, constrainTo),
                           names = names(modelTools$start))
        varPredictors <- modelTools$varPredictors(theta)
        eta <- offset + modelTools$predictor(varPredictors)
        if (!missing(eliminate)) eta <- eta + alpha[eliminate]
        if (any(!is.finite(eta))) {
            stop("Values of 'start' and 'constrain' produce non-finite ",
                 "predictor values")
        }
        mu <- family$linkinv(eta)
        dev[1] <- sum(family$dev.resids(y, mu, weights))
        if (trace)
            prattle("Initial Deviance = ", dev, "\n", sep = "")

        if (status == "not.converged") {
            X <-  modelTools$localDesignFunction(theta, varPredictors)
            X <- X[, !isConstrained, drop = FALSE]
            #X <- Matrix(X)
            #classX <- class(X)
            if (!missing(eliminate)){
                nelim <- nlevels(eliminate)
                elim <- seq_len(nelim) + 1
                non.elim <- seq.int(elim + 1, length(theta))
                Tvec <- 1 + rep.int(ridge, nelim)
                grp.size <- tabulate(eliminate)
                grp.end <- cumsum(grp.size)
            }
            ridge <- 1 + ridge
            for (iter in seq(iterMax)) {
                if (any(is.infinite(X))){
                #if (any(is.infinite(X@x))){
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
                z <- (abs(dmu) >= eps) * (y - mu)/dmu
                vmu <- family$variance(mu)
                w <- weights * (abs(dmu) >= eps) * dmu * dmu/vmu
                if (any(!is.finite(w))) {
                    status <- "w.not.finite"
                    break
                }

                wSqrt <- sqrt(w)
                W.X <- wSqrt * X
                w.z <- wSqrt * z
                score <- drop(crossprod(W.X, w.z))
                diagInfo <- colSums(W.X^2)
                Xscales <- sqrt(diagInfo)
                Xscales[Xscales < 1e-3] <- 1e-3 ## to allow for zeros
                W.X.scaled <- t(W.X)/Xscales
                if (!missing(eliminate)) {
                    elim.diagInfo <- grp.sum(w, grp.end)
                    diagInfo <- c(elim.diagInfo, diagInfo)
                    score <- c(grp.sum(wSqrt * w.z, grp.end), score)
                }
                if (all(diagInfo < 1e-20) ||
                    all(abs(score) <
                        tolerance * sqrt(tolerance + diagInfo))) {
                    status <- "converged"
                    break
                }
                znorm <- sqrt(sum(w.z^2))
                if (!missing(eliminate)){
                    w.z <- w.z/znorm
                    W.Z <- rbind2(w.z, W.X.scaled)
                    elimXscales <- sqrt(elim.diagInfo)
                    Umat <- rowsum(as.matrix(wSqrt*t(W.Z)), eliminate)/elimXscales
                    Wmat <- tcrossprod(W.Z)
                    diag(Wmat) <- ridge
                    ZWZinv <- solve1(Wmat, Tvec, Umat, elim)
                    alphaChange <- -ZWZinv[elim]/ZWZinv[1] * znorm/elimXscales
                    thetaChange <- -ZWZinv[non.elim]/ZWZinv[1] * znorm/Xscales
                } else {
                    XWX <- tcrossprod(W.X.scaled)
                    diag(XWX) <- ridge
                    theChange <- solve(XWX, score/(znorm * Xscales)) * znorm/Xscales
                }
                dev[2] <- dev[1]
                j <- 1
                while (dev[1] >= dev[2] && j < 11) {
                    nextTheta <- replace(theta, !isConstrained,
                                         theta[!isConstrained] + thetaChange)
                    varPredictors <- modelTools$varPredictors(nextTheta)
                    eta <- offset + modelTools$predictor(varPredictors)
                    if (!missing(eliminate)) {
                        nextAlpha <- alpha + alphaChange
                        eta <- eta + nextAlpha[eliminate]
                    }
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
                X <- modelTools$localDesignFunction(theta, varPredictors)
                X <- X[, !isConstrained, drop = FALSE]
                #X <- as(X, classX)
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
            if (attempt > 5)#|| !any(unspecifiedNonlin))
                return()
            else if (verbose)
                message("Restarting")
        }
    }
    theta[constrain] <- NA
    if (!missing(eliminate)) {
        ## sweeps needed to get the rank right
        subtracted <- rowsum(t(W.X.scaled), eliminate)/grp.size
        subtracted[,1] <- 0
        W.X.scaled <- t(W.X.scaled) - subtracted[eliminate,]
        XWX <- crossprod(W.X.scaled)
    }
    else if (is.null(XWX)) XWX <- tcrossprod(W.X.scaled)
    Svd <- svd(XWX, nu = 0, nv = 0)
    sv.tolerance <-  100 * .Machine$double.eps
    theRank <- sum(Svd$d > max(sv.tolerance * Svd$d[1], 0)) + nelim
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs),
                                            mu, weights, dev[1])
                                 + 2 * theRank)
    fit <- list(coefficients = c(alpha, theta), constrain = constrain,
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
        varPredictors <- modelTools$varPredictors(theta)
        fit$termPredictors <- modelTools$predictor(varPredictors, term = TRUE)
    }
    fit
}
