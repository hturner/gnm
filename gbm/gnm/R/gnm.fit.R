gnm.fit <- function(modelTools, y, constrain, family = poisson(),
                     weights = rep.int(1, length(y)),
                     offset = rep.int(0, length(y)), nObs = length(y),
                     start = NULL, control = gnm.control(...), x, vcov) {
    conv <- FALSE
    repeat{
        if (is.null(start)) {
            theta <- modelTools$start()
            theta[constrain] <- 0
            for (iter in seq(control$startit)) {
                for (relevant in split(seq(length(theta)),
                                       modelTools$factorAssign)) {
                    parameterList <- unname(split(theta,
                                                  modelTools$factorAssign))
                    factorList <- modelTools$factorList(parameterList)
                    eta <- offset + modelTools$predictor(factorList)
                    X <- modelTools$localDesignFunction(parameterList,
                                                        factorList)
                    mu <- family$linkinv(eta)
                    dmu <- family$mu.eta(eta)
                    vmu <- family$variance(mu)
                    score <- crossprod(weights*(y - mu)*dmu/vmu,
                                       X[,relevant])
                    gradient <- crossprod(weights*dmu*dmu/vmu,
                                          (X[,relevant])^2)
                    theta[relevant] <- ifelse(!constrain[relevant],
                                              as.vector(theta[relevant] +
                                                 score/gradient), 0)
                }
            }
        }    
        else theta <- ifelse(!constrain, start, 0)
        for (iter in seq(control$maxit)) {
            parameterList <- unname(split(theta, modelTools$factorAssign))
            factorList <- modelTools$factorList(parameterList)
            eta <- offset + modelTools$predictor(factorList)
            X <- modelTools$localDesignFunction(parameterList, factorList)
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
            score <- as.vector(crossprod(weights*(y - mu)*dmu/vmu, X))
            if (all(abs(score) < control$epsilon*sqrt(diag(Info)))){
                conv <- TRUE
                break
            }
            theta <- ifelse(!constrain,
                        theta + drop(crossprod(VCOV, score)), 0)
        }
        if (!inherits(VCOV, "try-error") | !is.null(start)) break
        else cat("Bad parameterisation, restarting.\n")
    }
    if (!conv)
        warning("Fitting algorithm has either not converged or converged to\n",
                "a non-solution of the likelihood equations.\n",
                "Re-start gnm with coefficients of returned model.\n")
    theta[constrain] <- NA
    modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs), mu, weights,
                                            dev) + 2 * attr(VCOV, "rank"))
    fit <- list(coefficients = structure(theta, names =
                names(modelTools$factorAssign)), predictors = eta,
                fitted.values = mu, deviance = dev, aic = modelAIC,
                iter = iter,
                conv = conv, weights = diag(W), residuals = (y - mu)/dmu,
                df.residual = nObs - attr(VCOV, "rank"),
                rank = attr(VCOV, "rank"))
    if (x) fit$x <- structure(X, assign = modelTools$termAssign, dimnames =
                              list(NULL, names(modelTools$factorAssign)))
    if (vcov) {
        VCOV[constrain, constrain] <- 0
        fit$vcov <- structure(VCOV, dimnames =
                      c(rep.int(list(names(modelTools$factorAssign)), 2)))
    }
    fit    
}

