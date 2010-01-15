##  Functions to help compute group sums/means quickly
quickdiff <- function(x) {x[-1] - x[-length(x)]}

grp.sum <- function(x, grp.end){
    quickdiff(c(0, cumsum(as.numeric(x))[grp.end]))
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
              accelerate = 5 ## either numeric > 2, or FALSE
              )
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
    if (ordTRUE <- !identical(ord, xtfrm(eliminate))) {
        y <- as.numeric(y[ord])
        weights <- weights[ord]
        offset <- offset[ord]
        x <- x[ord,]
        eliminate <- eliminate[ord]
    }
    size <- tabulate(eliminate)
    end <- cumsum(size)
    nelim <- nlevels(eliminate)
    if (is.null(start)) { # use either y or etastart or mustart
        if (is.null(mustart) && is.null(etastart)) {
            elim.means <- grp.sum(y, end)/size
            os.by.level <- link(0.999 * elim.means + 0.001 * mean(y))
        } else {
            if (!is.null(etastart)) mustart <- linkinv(etastart)
            os.by.level <- link(grp.sum(mustart, end)/size)
        }
    } else start[seq(nelim)]
    os.vec <- offset + os.by.level[eliminate]
    eta.stored <- eta <- os.vec
    mu <- linkinv(eta)
    mu.eta <- linkder(eta)
    z <- eta + (y - mu) / mu.eta
    w <- weights * (mu.eta)^2/variance(mu)
    counter <- 0
    devold <- 0
#    if (!intercept) {
#        x <- cbind(1, x)  ## this typically results in fewer iterations
#    }
    if (intercept) x <- x[, -1, drop = FALSE]
    ## sweeps needed to get the rank right
    subtracted <- rowsum(x, eliminate)/size
    subtracted[,1] <- 0
    x <- x - subtracted[eliminate,]
    for (i in 1:control$maxit) {
        model <- lm.wfit(x, z, w, offset = os.vec)
        eta <- model$fitted
        mu <- linkinv(eta)
        mu.eta <- linkder(eta)
        res <- (y - mu) / mu.eta
        w <- weights * (mu.eta)^2/variance(mu)
        os.by.level <- os.by.level + grp.sum(w * res, end)/grp.sum(w, end)
        if (accelerate > 2) {
            ## store the current and previous two iterates
            if (i > 2) os.by.level.lag2 <- os.by.level.lag1
            if (i > 1) os.by.level.lag1 <- os.by.level.kept
            os.by.level.kept <- os.by.level
            counter <- counter + 1
            ## if i is a multiple of accelerate, then make an Aitken
            ## prediction step (for each eliminated parameter)
            if (counter == accelerate) {
                diff1 <- os.by.level.lag1 - os.by.level.lag2
                diff0 <- os.by.level - os.by.level.lag1
                if (all(abs(diff0 - diff1) > 1e-12)) {
                    os.by.level <- os.by.level.lag2 - diff1^2/(diff0 - diff1)
                } else accelerate <- 0 ## stop doing Aitken steps if
                ## differences are too small
                counter <- 0
            }
        }
        old.os.vec <- os.vec
        os.vec <- as.vector(offset + os.by.level[eliminate])
        eta <- eta - old.os.vec + os.vec
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
        z <- eta + (y - mu) / mu.eta
        w <- weights * (mu.eta)^2/variance(mu)
    }
    converged <- !(i == control$maxit)
    if (!converged) warning(paste("The convergence criterion was not met after",
                                  control$maxit, "iterations."))
    rank <- (model$rank + nelim)
    aic.model <- aic(y, sum(weights > 0), mu, weights, dev) + 2 * rank
    na.to.zero <- function(vec) ifelse(is.na(vec), 0, vec)
    coefs.elim <- os.by.level - subtracted %*% na.to.zero(model$coef)
    if (ordTRUE) {
        reorder <- order(ord)
        y <- y[reorder]
        mu <- mu[reorder]
        eta <- eta[reorder]
        weights <- weights[reorder]
        offset <- offset[reorder]
    }
    mu.eta <- linkder(eta)
    z <- eta + (y - mu) / mu.eta
    w <- weights * (mu.eta)^2/variance(mu)
    list(coefficients = c(coefs.elim, model$coef),
         residuals = (y - mu) / mu.eta,
         fitted.values = mu,
         rank = rank,
         family = family,
         linear.predictors = eta,
         deviance = dev,
         aic = aic.model,
         null.deviance = sum(dev.resids(y, linkinv(offset), weights)),
         iter = i,
         weights = weights * (mu.eta)^2/variance(mu),
         prior.weights = weights,
         df.residual = nobs - sum(weights == 0) - rank,
         df.null = nobs - sum(weights == 0),
         y = y,
         converged = converged)
    ##  NB: some components of the result of glm.fit are missing from this list
}


