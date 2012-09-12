gnmFit <-
    function (modelTools, y,
              constrain = numeric(0), # index of non-eliminated parameters
              constrainTo = numeric(length(constrain)),
              eliminate = NULL, # now a factor
              family = poisson(),
              weights = rep.int(1, length(y)),
              offset = rep.int(0, length(y)),
              nobs = length(y),
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
    names(y) <- NULL
    eps <- 100*.Machine$double.eps
    attempt <- 1
    if (verbose)
        width <- as.numeric(options("width"))
    nTheta <- length(modelTools$start)
    nelim <- nlevels(eliminate)
    non.elim <- seq.int(nelim + 1, length(start))

    ## add constraints specified by modelTools and glm
    tmpTheta <- as.double(rep(NA, nTheta))
    varPredictors <- modelTools$varPredictors(tmpTheta)
    X <- modelTools$localDesignFunction(tmpTheta, varPredictors)
    isLinear <- unname(!is.na(colSums(X)))
    if (any(isLinear)) {
        tmpTheta[isLinear] <-
            suppressWarnings(glm.fit.e(X[, isLinear, drop = FALSE], y,
                                       family = family,
                                       intercept = FALSE,
                                       eliminate = if (nelim) eliminate
                                       else NULL,
                                       coefonly = TRUE,
                                       control = glm.control(maxit = 1)))

        extraLin <- which(isLinear & is.na(tmpTheta))
    } else extraLin <- numeric()

    extra <- setdiff(c(modelTools$constrain, extraLin), constrain)
    ind <- order(c(constrain, extra))
    constrain <- c(constrain, extra)[ind]
    constrainTo <- c(constrainTo, numeric(length(extra)))[ind]
    isConstrained <- is.element(seq(nTheta), constrain)
    repeat {
        status <- "not.converged"
        unspecifiedNonlin <- FALSE
        dev <- numeric(2)
        if (nelim)  {
            elim <- seq.int(nelim)
            alpha <- start[elim]
        }
        else {
            eliminate <- 1
            alpha <- 0
        }
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
            ## update any unspecified linear parameters
            unspecified <- unname(is.na(theta))
            unspecifiedLin <- unspecified & isLinear
            unspecifiedNonlin <- unspecified & !isLinear
            if (!is.null(mustart))
                etastart <- family$linkfun(mustart)
            if (any(unspecifiedNonlin) && is.null(etastart)){
                theta[unspecifiedNonlin] <- gnmStart(sum(unspecifiedNonlin))
            }
            if (any(unspecifiedLin) || initElim) {
                ## offset nonLin terms (currently NA if using etastart)
                ## plus offset contribution of any specified lin par
                if (!is.null(etastart)) z <-  family$linkinv(etastart)
                else z <- y
                varPredictors <- modelTools$varPredictors(theta)
                tmpOffset <- modelTools$predictor(varPredictors, term = TRUE)
                tmpOffset <- rowSums(naToZero(tmpOffset))
                tmpOffset <- offset + alpha[eliminate] + tmpOffset
                ## starting values for elim ignored here
                tmpTheta <- suppressWarnings({
                    glm.fit.e(X[, unspecifiedLin, drop = FALSE],
                              z,
                              weights = weights,
                              etastart = etastart,
                              offset = tmpOffset,
                              family = family,
                              intercept = FALSE,
                              eliminate = if (nelim) eliminate else NULL,
                              coefonly = TRUE)})
                ## if no starting values for elim, use result of above
                if (initElim) alpha <- unname(attr(tmpTheta, "eliminated"))
                theta[unspecifiedLin] <- naToZero(tmpTheta)
            }
            if (any(unspecifiedNonlin) && !is.null(etastart)){
                ## offset linear terms
                ## plus contribution of specified nonlin terms
                varPredictors <- modelTools$varPredictors(theta)
                tmpOffset <- modelTools$predictor(varPredictors, term = TRUE)
                tmpOffset <- rowSums(naToZero(tmpOffset))
                tmpOffset  <- offset + alpha[eliminate] + tmpOffset
                if (any(isLinear) && isTRUE(all.equal(unname(etastart), tmpOffset))) {
                    etastart <- mustart <- NULL
                    eval(family$initialize)
                    etastart <- family$linkfun(mustart)
                }
                tmpOffset <- offset + alpha[eliminate]
                rss <- function(par) {
                    theta[unspecifiedNonlin] <<- par
                    varPredictors <<- modelTools$varPredictors(theta)
                    eta <<- tmpOffset + modelTools$predictor(varPredictors)
                    sum((etastart - eta)^2)
                }
                gr.rss <- function(par) {
                    X <- modelTools$localDesignFunction(theta, varPredictors)
                    -2 * t(X[, unspecifiedNonlin]) %*% ((etastart - eta))
                }
                theta[unspecifiedNonlin] <- optim(gnmStart(sum(unspecifiedNonlin)),
                                                  rss, gr.rss,
                                                  method = c("L-BFGS-B"),
                                                  control = list(maxit = iterStart),
                                                  lower = -10, upper = 10)$par
            }
            varPredictors <- modelTools$varPredictors(theta)
            tmpOffset <- offset + alpha[eliminate]
            eta <- tmpOffset + modelTools$predictor(varPredictors)
            mu <- family$linkinv(eta)
            dev[1] <- sum(family$dev.resids(y, mu, weights))
            if (trace)
                prattle("Initial Deviance = ",
                        format(dev[1], nsmall = 6), "\n", sep = "")
            niter <- iterStart * (any(unspecifiedNonlin) && is.null(etastart))
            for (iter in seq(length = niter)) {
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
                    Xi <- modelTools$localDesignFunction(theta,
                                                         varPredictors, i)
                    wXi <- weights * (abs(dmu) >= eps) * dmu * dmu/vmu * Xi
                    score <- sum((abs(y - mu) >= eps) * (y - mu)/dmu * wXi)
                    gradient <- sum(wXi * Xi)
                    otheta <- theta[i]
                    theta[i] <- as.vector(theta[i] + score/gradient)
                    if (!is.finite(theta[i])) {
                        status <- "bad.param"
                        break
                    }
                    varPredictors <- modelTools$varPredictors(theta)
                    eta <- tmpOffset + modelTools$predictor(varPredictors)
                    mu <- family$linkinv(eta)
                    ## poor man's step-halving
                    #if (sum(family$dev.resids(y, mu, weights)) > dev[1]){
                    #    theta[i] <- otheta
                    #    varPredictors <- modelTools$varPredictors(theta)
                    #    eta <- tmpOffset + modelTools$predictor(varPredictors)
                    #    mu <- family$linkinv(eta)
                    #}
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
                        alpha <- unname(attr(tmpTheta, "eliminated"))
                        tmpOffset <- offset + alpha[eliminate]
                    }
                    theta[which] <- tmpTheta
                    varPredictors <- modelTools$varPredictors(theta)
                    eta <- tmpOffset + modelTools$predictor(varPredictors)
                    mu <- family$linkinv(eta)
                }
                dev[1] <- sum(family$dev.resids(y, mu, weights))
                if (trace)
                    prattle("Start-up iteration ", iter, ". Deviance = ",
                            format(dev[1], nsmall = 6), "\n", sep = "")
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
                prattle("Initial Deviance = ", format(dev[1], nsmall = 6), "\n", sep = "")
        }
        if (status == "not.converged") {
            X <-  modelTools$localDesignFunction(theta, varPredictors)
            X <- X[, !isConstrained, drop = FALSE]
            np <- ncol(X) + 1
            ZWZ <- array(dim = c(np, np))
            I1 <- numeric(np)
            I1[1] <- 1
            if (nelim) Umat <- array(dim = c(nelim, np))
            #X <- Matrix(X)
            #classX <- class(X)
            if (nelim){
                grp.size <- tabulate(eliminate)
                grp.end <- cumsum(grp.size)
            }
            tmpAlpha <- 0
            ridge <- 1 + ridge
            for (iter in seq(length.out = iterMax + 1)) {
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
                w <- sqrt(weights * (abs(dmu) >= eps) * dmu * dmu/vmu)

                X <- w * X
                z <- w * z
                ZWZ[-1,-1] <- crossprod(X)
                score <- ZWZ[1,-1] <- ZWZ[-1,1] <- crossprod(z, X)
                ZWZ[1,1] <- sum(z * z)
                diagInfo <- diag(ZWZ)
                ## only check for non-eliminated coefficients
                if (any(!is.finite(diagInfo))) {
                    status <- "fail"
                    break
                }
                if (all(diagInfo < 1e-20) ||
                    all(abs(score) <
                        tolerance * sqrt(tolerance + diagInfo[-1]))) {
                    status <- "converged"
                    break
                }
                Zscales <- sqrt(diagInfo)
                Zscales[Zscales < 1e-3] <- 1e-3 ## to allow for zeros
                if (iter > iterMax) break
                if (nelim){
                    elimXscales <- grp.sum(w * w, grp.end)
                    elimXscales <- sqrt(elimXscales * ridge)
                    Umat[,1] <- rowsum.default(w * z, eliminate, reorder = FALSE)
                    Umat[,-1] <- rowsum.default(w * X, eliminate,
                                                reorder = FALSE)
                    Umat <- Umat/(elimXscales %o% Zscales)
                    ZWZ <- ZWZ/(Zscales %o% Zscales)
                    diag(ZWZ) <- ridge
                    z <- solve(ZWZ - crossprod(Umat), I1,
                                         tol = .Machine$double.eps)
                    thetaChange <- -z[-1]/z[1] * Zscales[1]/Zscales[-1]
                    alphaChange <- c(Umat %*% (z * sqrt(ridge)))/z[1] *
                                    Zscales[1]/elimXscales
                } else {
                    ZWZ <- ZWZ/(Zscales %o% Zscales)
                    diag(ZWZ) <- ridge
                    z <- solve(ZWZ, I1, tol = .Machine$double.eps)/Zscales
                    thetaChange <- -z[-1]/z[1]
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
                    mu <- family$linkinv(eta)
                    dev[1] <- sum(family$dev.resids(y, mu, weights))
                    if (!is.finite(dev[1])) {
                        status <- "no.deviance"
                        break
                    }
                    scale <- scale*2
                    j <- j + 1
                }
                if (status == "no.deviance") break
                if (trace){
                    prattle("Iteration ", iter, ". Deviance = ",
                            format(dev[1], nsmall = 6),
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
            if (any(!is.finite(eta)))
                 status <- "eta.not.finite"
             if (any(!is.finite(w)))
                 status <- "w.not.finite"
             if (any(is.infinite(X)))
                 status <- "X.not.finite"

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
                               no.deviance = "Deviance is not finite"))
            attempt <- attempt + 1
            if (attempt > 5 || (any(is.na(start)) && !any(unspecifiedNonlin)))
                return()
            else if (verbose)
                message("Restarting")
            X <- modelTools$localDesignFunction(theta, varPredictors)
            ridge <- ridge - 1
        }
    }
    theta[constrain] <- NA
    X <- modelTools$localDesignFunction(theta, varPredictors)
    X <- X[, !isConstrained, drop = FALSE]
    if (nelim) {
        ## sweeps needed to get the rank right
        subtracted <- rowsum.default(X, eliminate, reorder = FALSE)/grp.size
        if (modelTools$termAssign[1] == 0) subtracted[,1] <- 0
        theRank <- rankMatrix(X - subtracted[eliminate,]) + nelim
        names(alpha) <- paste("(eliminate)", elim, sep = "")
    }
    else theRank <- rankMatrix(X)
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nobs),
                                            mu, weights, dev[1])
                                 + 2 * theRank)
    fit <- list(coefficients = structure(theta, eliminated = alpha),
                constrain = constrain,
                constrainTo = constrainTo, residuals = z/w,
                fitted.values = mu,
                rank = theRank, family = family, predictors = eta,
                deviance = dev[1], aic = modelAIC,
                iter = iter - (iter != iterMax), weights = w * w,
                prior.weights = weights,
                df.residual = c(nobs - theRank),
                y = y)
    if (status == "not.converged") {
        warning("Fitting algorithm has either not converged or converged\n",
                "to a non-solution of the likelihood equations.\n",
                "Use exitInfo() for numerical details of last iteration.\n")
        fit$converged <- structure(FALSE, score = score, criterion =
                                   tolerance * sqrt(tolerance + diagInfo[-1]))
    }
    else
        fit$converged <- TRUE

    if (x) {
        X <- modelTools$localDesignFunction(theta, varPredictors)
        fit$x <- structure(X, assign = modelTools$termAssign)
    }
    if (termPredictors) {
        theta[is.na(theta)] <- 0
        varPredictors <- modelTools$varPredictors(theta)
        fit$termPredictors <- modelTools$predictor(varPredictors, term = TRUE)
    }
    fit
}
