profile.gnm <- function (fitted, which = 1:p, alpha = 0.01, maxsteps = 10, del = zmax/5, 
    trace = FALSE, ...) 
{
    Pnames <- names(B0 <- coefficients(fitted))
    pv0 <- t(as.matrix(B0))
    p <- length(Pnames)
    rnk <- fitted$rank
    if (is.character(which)) 
        which <- match(which, Pnames)
    summ <- summary(fitted)
    std.err <- summ$coefficients[, "Std. Error"]
    mf <- update(fitted, method = "model.frame")
    n <- length(Y <- model.response(mf))
    resdf <- fitted$df.residual
    O <- model.offset(mf)
    if (!length(O)) 
        O <- rep(0, n)
    W <- model.weights(mf)
    if (length(W) == 0) 
        W <- rep(1, n)
    OriginalDeviance <- deviance(fitted)
    DispersionParameter <- summ$dispersion
    fam <- family(fitted)
    switch(fam$family, binomial = {
        if (!is.null(dim(Y))) {
            n <- n/2
            O <- O[1:n]
            Y <- Y[, 1]/(W <- drop(Y %*% c(1, 1)))
        }
        zmax <- sqrt(qchisq(1 - alpha/2, rnk))
        profName <- "z"
    }, poisson = , "Negative Binomial" = {
        zmax <- sqrt(qchisq(1 - alpha/2, rnk))
        profName <- "z"
    }, gaussian = , quasi = , inverse.gaussian = , quasibinomial = , 
        quasipoisson = , {
            zmax <- sqrt(rnk * qf(1 - alpha/2, rnk, resdf))
            profName <- "tau"
        })
    origConstrain <- fitted$constrain
    prof <-  as.list(rep(NA, length(which)))
    names(prof) <- Pnames[which]
    which <- which[!is.na(std.err)[which]]
    keep.del <- del
    for (i in which) {
        zi <- 0
        pvi <- pv0
        pi <- Pnames[i]
        for (sgn in c(-1, 1)) {
            if (trace) 
                cat("\nParameter:", pi, c("down", "up")[(sgn + 
                  1)/2 + 1], "\n")
            step <- 0
            z <- 0
            init <- coef(fitted)
            del <- keep.del
            keep.pvi <- pvi
            keep.zi <- zi
            check <- 0
            while ((step <- step + 1) < maxsteps && abs(z) < 
                zmax) {
                bi <- B0[i] + sgn * step * del * std.err[i]
                fm <- try(update(fitted, constrain = rbind(origConstrain,
                                         data.frame(constrain = i, value = bi)),
                                 trace = FALSE, verbose = FALSE,
                                 start = init),
                          silent = TRUE)
                if (is.null(fm)) {
                    message("Could not complete profile for", pi, "\n")
                    break
                }
                init <- coef(fm)
                ri <- pv0
                ri[, names(coef(fm))] <- coef(fm)
                ri[, pi] <- bi
                pvi <- rbind(pvi, ri)
                zz <- (fm$deviance - OriginalDeviance)/DispersionParameter
                if (zz > -0.001) 
                  zz <- max(zz, 0)
                else stop("profiling has found a better solution, so original fit had not converged")
                z <- sgn * sqrt(zz)
                zi <- c(zi, z)
                print(data.frame(step = step, val = bi, deviance = fm$deviance,
                                 zstat = z))
            }
        }
        si <- order(zi)
        prof[[pi]] <- structure(data.frame(zi[si]), names = profName)
        prof[[pi]]$par.vals <- pvi[si, ]
    }
    val <- structure(prof, original.fit = fitted, summary = summ)
    class(val) <- c("profile.gnm", "profile.glm", "profile")
    val
}
