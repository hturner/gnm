rowsum.unique <- function (x, group, ugroup,...)
{
    rval <- .Call("Rrowsum_matrix", x, NCOL(x), group, ugroup,
        FALSE, PACKAGE = "base")
    dimnames(rval) <- list(ugroup, colnames(x))
    rval
}
"glm.fit.e" <-
##  This fits a glm with eliminated factor, and should be much quicker
##  than glm.fit when the number of levels of the eliminated factor is large.
##  It is designed to work with the slightly modified glm function, 'glm.e'.
##
##  The 'eliminate' argument, if not NULL, must be a factor.
##
##  No account is taken of NAs -- will that be a problem, or have they gone by
##  the time glm.fit.e gets called?  None of the other extensive checks done by
##  glm.fit are done here either.
##
##  Proof of concept only, no guarantee of fitness for any purpose.
##
##  David Firth, September 2009.
##
    function (x, y,
              weights = rep(1, NROW(y)),
              start = NULL,
              etastart = NULL,
              mustart = NULL,
              offset = rep(0, NROW(y)),
              family = gaussian(),
              control = glm.control(), ## only for compatibility with glm.fit
              intercept = TRUE, ## only for compatibility with glm.fit
              eliminate = NULL,  ## alternatively a factor
              ridge = 1e-8)
{
    if (is.null(eliminate)) { ## just revert to glm.fit
        return(glm.fit(x, y, weights = weights, start = start,
                       etastart = etastart, mustart = mustart,
                       offset = offset, family = family,
                       control = control,
                       intercept = intercept))
    }
##  The rest handles the case of an eliminated factor
    nobs <- NROW(y)
    non.elim <- ncol(x)
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    link <- family$linkfun
    linkinv <- family$linkinv
    linkder <- family$mu.eta
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    ## sort data to help compute group means quickly
    ord <- order(xtfrm(eliminate))
    if (ordTRUE <- !identical(ord, seq(eliminate))) {
        y <- as.numeric(y[ord])
        weights <- weights[ord]
        offset <- offset[ord]
        if (non.elim) x <- x[ord,]
        eliminate <- eliminate[ord]
    }
    size <- tabulate(eliminate)
    end <- cumsum(size)
    nelim <- nlevels(eliminate)
    elim <- seq.int(nelim)
    if (is.null(start)) { # use either y or etastart or mustart
        if (is.null(mustart) && is.null(etastart)) {
            elim.means <- grp.sum(y, end)/size
            os.by.level <- link(0.999 * elim.means + 0.001 * mean(y))
        } else {
            if (!is.null(etastart)) mustart <- linkinv(etastart)
            os.by.level <- link(grp.sum(mustart, end)/size)
        }
    } else os.by.level <- start[elim]
    os.vec <- offset + os.by.level[eliminate]
    eta.stored <- eta <- os.vec
    mu <- linkinv(eta)
    mu.eta <- linkder(eta)
    z <- eta + (y - mu) / mu.eta
    w <- weights * (mu.eta)^2/variance(mu)
    counter <- 0
    devold <- 0
    if (intercept) x <- x[, -1, drop = FALSE]
    if (non.elim) {
        ## sweeps needed to get the rank right
        subtracted <- rowsum.unique(x, eliminate, elim)/size
        subtracted[,1] <- 0
        x <- x - subtracted[eliminate,]
        ## initial fit to drop aliased columns
        model <- lm.wfit(x, z, w, offset = os.vec)
        eta <- model$fitted
        mu <- linkinv(eta)
        mu.eta <- linkder(eta)
        z <- eta + (y - mu) / mu.eta
        w <- weights * (mu.eta)^2/variance(mu)
        full.theta <- model$coefficients
        est <- !is.na(full.theta)
        x <- x[, est]
        theta <- full.theta[est]
    }
    Z <- cbind(z, x)
    I1 <- numeric(ncol(Z))
    I1[1] <- 1
    for (i in 1:control$maxit) {
        ## try without scaling etc - already of full rank
        W.Z <- w * Z
        Tvec <- grp.sum(w, end)
        Umat <- rowsum.unique(W.Z, eliminate, elim)
        Wmat <- crossprod(Z, W.Z)
        diag(Wmat) <- diag(Wmat) + ridge
        Ti.U <- Umat/Tvec
        Qi <- solve(Wmat - crossprod(Umat, Ti.U), I1)
        theta <- -Qi[-1]/Qi[1]
        os.by.level <- (Ti.U %*% Qi)/Qi[1]

        if (non.elim) eta <- drop(x %*% theta + offset + os.by.level[eliminate])
        else eta <- offset + os.by.level[eliminate]
        mu <- linkinv(eta)
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace)
            cat("Deviance =", dev, "Iterations -", i,
                "\n")
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
            conv <- TRUE
            coef <- start
            break
        }
        devold <- dev
        mu.eta <- linkder(eta)
        Z[,1] <- eta + (y - mu) / mu.eta
        w <- weights * (mu.eta)^2/variance(mu)
    }
    converged <- !(i == control$maxit)
    if (!converged) warning(paste("The convergence criterion was not met after",
                                  control$maxit, "iterations."))
    if (non.elim) {
        rank <- (model$rank + nelim)
        full.theta[est] <- theta
        os.by.level <- os.by.level - subtracted %*% naToZero(full.theta)
    }
    else {
        rank <- nelim
        full.theta <- numeric(0)
    }
    names(os.by.level) <- paste("(eliminate)", elim, sep = "")
    aic.model <- aic(y, sum(weights > 0), mu, weights, dev) + 2 * rank
    if (ordTRUE) {
        reorder <- order(ord)
        y <- y[reorder]
        mu <- mu[reorder]
        eta <- eta[reorder]
        weights <- weights[reorder]
        offset <- offset[reorder]
    }
    mu.eta <- linkder(eta)
    list(coefficients = c(os.by.level, full.theta),
         residuals = (y - mu) / mu.eta,
         fitted.values = mu,
         rank = rank,
         family = family,
         predictors = eta,
         deviance = dev,
         aic = aic.model,
         iter = i,
         weights = weights * (mu.eta)^2/variance(mu),
         prior.weights = weights,
         df.residual = nobs - sum(weights == 0) - rank,
         y = y,
         converged = converged)
    ##  NB: some components of the result of glm.fit are missing from this list
}


