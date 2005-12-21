"quick.glm.fit" <-
##  A wrapper for glm.fit, which is much faster when a large number
##  of parameters can be eliminated, but which typically (if nIter is small)
##  stops before convergence.  Useful for getting gnm starting values.
##
##  The eliminate argument is assumed numeric (indices of columns in X);
##  or logical (indicating for each column of x whether to eliminate).  No
##  check is done here on "eliminability" of the specified columns.
##
##  The non-eliminated columns are assumed not to include the intercept (ie
##  no column of ones).
##
##  When eliminate is used, only the "coefficients" component is returned.
##  (for reasons of speed/laziness).  This is fine for gnm purposes, but if
##  quick.glm.fit is made into a `method' for glm() fits then the result
##  needs to have various other components added.
##
##  No account is taken of NAs -- will that be a problem, or have they gone by
##  the time glm.fit gets called?
##
    function (x, y,
              weights = rep(1, nobs),
              offset = rep(0, nobs),
              family = gaussian(),
              eliminate = NULL,
              nIter = 2,
              verbose = FALSE)
{
    if (is.null(eliminate)) return(glm.fit(x, y,
                                           weights = weights,
                                           offset = offset,
                                           family = family))
##  The rest handles the case of eliminated columns in X
    x <- as.matrix(x)
    if (is.logical(eliminate)){
        if (length(eliminate) != ncol(x))
            stop("eliminate argument has wrong length")
        eliminate <- which(eliminate)
    }
    xElim <- x[ , eliminate, drop = FALSE]
    xNotElim <- cbind(1, x[ , -eliminate, drop = FALSE])
    os.by.level <- rep(0, length(eliminate))
    model <- glm.fit(xNotElim, y,
                     weights = weights,
                     offset = offset,
                     family = family,
                     control = glm.control(maxit = 1))
    for (i in 1:nIter) {
        if (verbose) cat("quick.glm.fit iteration", i,
                         "deviance =", deviance(model), "\n")
        eta <- model$linear.predictors
        mu <- fitted(model)
        dmu <- family$mu.eta(eta)
        z <- model$residuals
        w <- model$weights
        w <- xElim * w
        wz <- w * z
        os.by.level <- os.by.level + colSums(wz)/colSums(w) + coef(model)[1]
        os.vec <- colSums(os.by.level * t(xElim))
        ow <- options("warn")
        options(warn = -1)
        model <- glm.fit(xNotElim, y,
                         weights = weights,
                         offset = offset + os.vec,
                         etastart = eta,
                         family = family,
                         control = glm.control(maxit = 2))
        options(ow)
    }
    coefs <- rep(NA, ncol(X))
    names(coefs) <- colnames(X)
    coefs[eliminate] <- os.by.level
    coefs[-eliminate] <- coef(model)[-1] ## intercept is dropped
    list(coefficients = coefs)
}
