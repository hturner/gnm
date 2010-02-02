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
              start = rep.int(NA, length(modelTools$start) + nlevels(eliminate)),
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
    nTheta <- length(modelTools$start)
    nelim <- nlevels(eliminate)
    if (nelim)  {
        elim <- seq.int(nelim)
        alpha <- start[elim]
        names(alpha) <- paste("(eliminate)", elim, sep = "")
    }
    else {
        eliminate <- 1
        alpha <- 0
    }
    non.elim <- seq.int(nelim + 1, length(start))
    isConstrained <- is.element(seq(nTheta), constrain)
    XWX <- NULL
    repeat {
        status <- "not.converged"
        unspecifiedNonlin <- FALSE
        if (any(is.na(start))) {
            if (verbose == TRUE)
                prattle("Initialising", "\n", sep = "")
            ## only use start for elim par if all specified
            initElim <- any(is.na(alpha))
            if (initElim) alpha[] <- numeric(nelim)
            theta <- start[non.elim]
            theta[is.na(theta)] <- modelTools$start[is.na(theta)]
            names(theta) <- names(modelTools$start)
            theta[constrain] <- constrainTo
            tmpTheta <- as.double(rep(NA, nTheta))
            varPredictors <- modelTools$varPredictors(tmpTheta)
            X <- modelTools$localDesignFunction(tmpTheta, varPredictors)
            ## update any unspecified linear parameters
            isLinear <- !is.na(colSums(X))
            unspecified <- is.na(theta)
            unspecifiedLin <- unspecified & isLinear
            unspecifiedNonlin <- unspecified & !isLinear
            if (any(unspecifiedLin) || !initElim) {
                ## offset any nonLin terms fully specified by start/modelTools$start/constrain
                ## plus offset contribution of any specified lin par
                theta[unspecifiedLin] <- 0
                varPredictors <- modelTools$varPredictors(theta)
                tmpOffset <- modelTools$predictor(varPredictors, term = TRUE)
                tmpOffset <- rowSums(naToZero(tmpOffset))
                tmpOffset <- offset + alpha[eliminate] + tmpOffset
                ## assume either elim all specified or all not specified
                tmpTheta <- suppressWarnings(glm.fit.e(X[, unspecifiedLin], y,
                                                       weights = weights,
                                                       offset = tmpOffset,
                                                       family = family,
                                                       intercept = FALSE,
                                                       eliminate =
                                                       if (initElim) eliminate
                                                       else NULL,
                                                       control = glm.control(maxit = 15))$coefficients)
                if (initElim) {
                    alpha[] <- tmpTheta[elim]
                    theta[unspecifiedLin] <- tmpTheta[-elim]
                }
                else theta[unspecifiedLin] <- tmpTheta
                if (sum(is.na(theta[isLinear])) > length(constrain)) {
                    extra <- setdiff(which(is.na(theta[isLinear])), constrain)
                    isConstrained[extra] <- TRUE
                    ind <- order(c(constrain, extra))
                    constrain <- c(constrain, extra)[ind]
                    constrainTo <- c(constrainTo, numeric(length(extra)))[ind]
                }
                theta <- naToZero(theta)
            }
            if (any(unspecifiedNonlin)){
                theta[unspecifiedNonlin] <- gnmStart(sum(unspecifiedNonlin))
            }
            if (!is.null(mustart))
                etastart <- family$linkfun(mustart)
            if (!is.null(etastart)){
                varPredictors <- modelTools$varPredictors(theta)
                tmpOffset <- modelTools$predictor(varPredictors, term = TRUE)
                tmpOffset <- rowSums(naToZero(tmpOffset))
                tmpOffset  <- offset + alpha[eliminate] + tmpOffset
                if (any(isLinear) && isTRUE(all.equal(unname(etastart), tmpOffset))) {
                    etastart <- start <- mustart <- NULL
                    eval(family$initialize)
                    etastart <- family$linkfun(mustart)
                }
                tmpTheta <- numeric(nTheta)
                rss <- function(theta) {
                    tmpTheta[unspecifiedNonlin] <<- theta
                    varPredictors <<- modelTools$varPredictors(tmpTheta)
                    eta <<- tmpOffset + modelTools$predictor(varPredictors)
                    sum((etastart - eta)^2)
                }
                gr.rss <- function(theta) {
                    X <- modelTools$localDesignFunction(theta, varPredictors)
                    -2 * t(X[, unspecifiedNonlin]) %*% ((etastart - eta))
                }
                theta[unspecifiedNonlin] <- optim(theta[unspecifiedNonlin],
                                                  rss, gr.rss,
                                                  method = c("L-BFGS-B"),
                                                  control = list(maxit = iterStart),
                                                  lower = -10, upper = 10)$par
                iterStart <- 0
            }
            varPredictors <- modelTools$varPredictors(theta)
            tmpOffset <- offset + alpha[eliminate]
            eta <- tmpOffset + modelTools$predictor(varPredictors)
            mu <- family$linkinv(eta)
            dev[1] <- sum(family$dev.resids(y, mu, weights))
            if (trace)
                prattle("Initial Deviance = ", dev[1], "\n", sep = "")
            for (iter in seq(length = iterStart * any(unspecifiedNonlin))) {
                if (verbose) {
                    if (iter == 1)
                        prattle("Running start-up iterations", "\n"[trace],
                                sep = "")
                    if ((iter + 25)%%width == (width - 1))
                        cat("\n")
                }
                for (i in rep.int(seq(nTheta)[unspecifiedNonlin], 2)) {
                    dmu <- family$mu.eta(eta)
                    vmu <- family$variance(mu)
                    w <- weights * (abs(dmu) >= eps) * dmu * dmu/vmu
                    Xi <- modelTools$localDesignFunction(theta,
                                                         varPredictors, i)
                    score <- crossprod((abs(y - mu) >= eps) * (y - mu)/dmu,
                                       w * Xi)
                    gradient <- crossprod(w, Xi^2)
                    theta[i] <- as.vector(theta[i] + score/gradient)
                    if (!is.finite(theta[i])) {
                        status <- "bad.param"
                        break
                    }
                    varPredictors <- modelTools$varPredictors(theta)
                    eta <- tmpOffset + modelTools$predictor(varPredictors)
                    mu <- family$linkinv(eta)
                }
                if (status == "not.converged" && any(isLinear)) {
                    if (iter == 1) {
                        which <- which(isLinear & !isConstrained)
                        if(!exists("X"))
                            X <- modelTools$localDesignFunction(theta,
                                                                varPredictors)
                    }
                    tmpTheta <- updateLinear(which, theta, y, mu, eta, offset,
                                             weights, family, modelTools, X,
                                             if(nelim) eliminate else NULL)
                    if (nelim){
                        alpha[] <- tmpTheta[elim]
                        theta[which] <- tmpTheta[-elim]
                        tmpOffset <- offset + alpha[eliminate]
                    }
                    else theta[which] <- tmpTheta
                    varPredictors <- modelTools$varPredictors(theta)
                    eta <- tmpOffset + modelTools$predictor(varPredictors)
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
            theta <- structure(replace(start[non.elim], constrain, constrainTo),
                               names = names(modelTools$start))
            varPredictors <- modelTools$varPredictors(theta)
            eta <- offset + alpha[eliminate] + modelTools$predictor(varPredictors)
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
            X <-  modelTools$localDesignFunction(theta, varPredictors)
            X <- X[, !isConstrained, drop = FALSE]
            #X <- Matrix(X)
            #classX <- class(X)
            if (nelim){
                Tvec <- 1 + rep.int(ridge, nelim)
                grp.size <- tabulate(eliminate)
                grp.end <- cumsum(grp.size)
            }
            tmpAlpha <- 0
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
                if (nelim) {
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
                if (nelim){
                    w.z <- w.z/znorm
                    W.Z <- rbind2(w.z, W.X.scaled)
                    elimXscales <- sqrt(elim.diagInfo)
                    Umat <- rowsum(as.matrix(wSqrt*t(W.Z)), eliminate)/elimXscales
                    Wmat <- tcrossprod(W.Z)
                    diag(Wmat) <- ridge
                    coef <- solve1(Wmat, Tvec, Umat, elim)
                    alphaChange <- coef[elim] * znorm/elimXscales
                    thetaChange <- coef[-elim] * znorm/Xscales
                } else {
                    XWX <- tcrossprod(W.X.scaled)
                    diag(XWX) <- ridge
                    thetaChange <- solve(XWX, score/(znorm * Xscales)) * znorm/Xscales
                }
                dev[2] <- dev[1]
                j <- scale <- 1
                while (dev[1] >= dev[2] && j < 11) {
                    if (nelim) tmpAlpha <- alpha + alphaChange/scale
                    tmpTheta <- replace(theta, !isConstrained,
                                        theta[!isConstrained] + thetaChange/scale)
                    varPredictors <- modelTools$varPredictors(tmpTheta)
                    eta <- offset + tmpAlpha[eliminate] +
                        modelTools$predictor(varPredictors)
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
                    scale <- scale*2
                    j <- j + 1
                }
                if (trace){
                    prattle("Iteration ", iter, ". Deviance = ", dev[1],
                            "\n", sep = "")
                }
                else if (verbose)
                    prattle(".")
                if (nelim) alpha <- tmpAlpha
                theta <- tmpTheta
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
    if (nelim) {
        ## sweeps needed to get the rank right
        subtracted <- rowsum(t(W.X.scaled), eliminate)/grp.size
        subtracted[,1] <- 0
        W.X.scaled <- t(W.X.scaled) - subtracted[eliminate,]
        XWX <- crossprod(W.X.scaled)
    }
    else {
        alpha <- numeric(0)
        if (is.null(XWX)) XWX <- tcrossprod(W.X.scaled)
    }
    diag(XWX) <- diag(XWX) - ridge + 1
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
