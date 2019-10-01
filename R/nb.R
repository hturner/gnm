#  A generic version of MASS::glm.nb.
#
#  Copyright (C) 1994-2014 W. N. Venables and B. D. Ripley
#  Copyright (C) 2019 Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# generic to work on fitted poisson/NB models inheriting from glm
import::from(stats, glm.control)
nb <- function(object, control = glm.control(...), ...){
    stopifnot(inherits(object, "glm"))
    stopifnot(object$family$family %in% c("poisson", "negative.binomial"))
    UseMethod()
}

# default uses update
import::from(stats, glm.control)
nb.default <- function(object, control = glm.control(...), ...){
    fit <- nb_fit(fit = object, fun = "update", control = control)
    fit$call <- match.call()
    fit$call[[1]] <- as.name("nb")
    class(fit) <- c("negbin", class(fit))
}

# overwrite summary and vcov method from MASS
summary.negbin <- function(object, dispersion = 1, correlation = FALSE, ...){
    if (is.null(dispersion)) dispersion <- 1
    summ1 <- NextMethod(dispersion = dispersion, correlation = correlation, ...)
    summ2 <- c(summ1,
               object[c("theta", "SE.theta", "twologlik", "th.warn")])
    class(summ2) <- c("summary.negbin", class(summ1))
    summ2
}

vcov.negbin <- function(object, dispersion = 1, ...){
    if (is.null(dispersion)) dispersion <- 1
    NextMethod(dispersion = dispersion, ...)
}

# generalized function to estimate theta 
import::from(stats, glm.control)
import::from(MASS, negative.binomial)
import::from(MASS, theta.ml)
nb.fit <- function(fit, # initial fit
                   fun, # fitting function
                   fun.args = NULL, # further arguments to function
                   control = glm.control(...), 
                   ...){ # other args to pass to fitting function
    
    # family should be Poisson or NB with specified theta
    stopifnot(fit$family$family %in% c("poisson", "negative.binomial"))
    fit.control <- list(maxit = control$maxit, 
                        epsilon = control$epsilon, 
                        trace = control$trace > 1)
    
    # estimate initial theta from initial fit
    y <- fit$y
    mu <- fit$fitted.values
    w <- fit$prior.weights
    th <- as.vector(theta.ml(y, mu, sum(w), w, 
                             limit = control$maxit, 
                             trace = control$trace > 2))
    if (control$trace) 
        message(gettextf("Initial value for 'theta': %f", signif(th)), 
                domain = NA)
    fam <- negative.binomial(theta = th, link = family$link)
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    Lm <- loglik(n, th, mu, y, w)
    Lm0 <- Lm + 2 * d1
    while ((iter <- iter + 1) <= control$maxit && 
           (abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon) {
        start <- fit$coefficients
        start[aliased] <- 0
        fit <- do.call(fun, c(list(family = fam, start = start, 
                                   control = fit.control), fun.args))
        t0 <- th
        mu <- fit$fitted.values
        th <- theta.ml(y, mu, sum(w), w, limit = control$maxit, 
                       trace = control$trace > 2)
        fam <- negative.binomial(theta = th, link = family$link)
        del <- t0 - th
        Lm0 <- Lm
        Lm <- loglik(n, th, mu, y, w)
        if (control$trace) {
            Ls <- loglik(n, th, y, y, w)
            Dev <- 2 * (Ls - Lm)
            message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f", 
                            iter, signif(th), signif(Dev)), domain = NA)
        }
    }
    if (!is.null(attr(th, "warn"))) 
        fit$th.warn <- attr(th, "warn")
    if (iter > control$maxit) {
        warning("alternation limit reached")
        fit$th.warn <- gettext("alternation limit reached")
    }
    # Cannot find example of case where extra iterations needed to compute 
    # null.deviance.
    fit$theta <- as.vector(th)
    fit$SE.theta <- attr(th, "SE")
    fit$twologlik <- as.vector(2 * Lm)
    fit$aic <- -fit$twologlik + 2 * fit$rank + 2
    fit
}
