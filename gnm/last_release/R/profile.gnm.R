profile.gnm <- function (fitted, which = ofInterest(fitted), alpha = 0.05,
                         maxsteps = 10, stepsize = NULL, trace = FALSE, ...)
{
    fittedCoef <- parameters(fitted)
    coefNames <- names(fittedCoef)
    p <- length(coefNames)
    if (is.null(which))
        which <- 1:p
    else if (is.numeric(which))
        which <- which
    else if (is.character(which))
        which <- match(which, coefNames)
    summ <- summary(fitted)
    sterr <- summ$coefficients[, "Std. Error"]
    fittedDev <- deviance(fitted)
    disp <- summ$dispersion
    ## use z cutoffs as in confint.profile.gnm
    zmax <- abs(qnorm(alpha/2))
    fittedConstrain <- fitted$constrain
    fittedConstrainTo <- fitted$constrainTo
    auto <- is.null(stepsize)
    if (!auto)
        stepsize[1:2] <- stepsize
    prof <-  as.list(rep(NA, length(which)))
    names(prof) <- coefNames[which]
    which <- which[!is.na(sterr)[which]]
    for (i in which) {
        par <- coefNames[i]
        prof[[par]] <- numeric(2 * maxsteps + 1)
        par.vals <- matrix(nrow = 2 * maxsteps + 1, ncol = p,
                           dimnames = list(NULL, coefNames))
        par.vals[maxsteps + 1,] <- fittedCoef
        asymptote <- c(FALSE, FALSE)
        if (auto) {
            ## set defaults
            sub <- 3 # no. of steps from MLE to zmax*se
            stepsize <- c(zmax/sub * sterr[i], zmax/sub * sterr[i])
            ## estimate quadratic in the region MLE +/- zmax*se
            margin <- zmax * sterr[i]
            updatedDev <- numeric(2)
            for (sgn in c(-1, 1)) {
                val <- fittedCoef[i] + sgn * margin
                updated <-
                    suppressWarnings(update(fitted, constrain =
                                            c(fittedConstrain, i),
                                            constrainTo =
                                            c(fittedConstrainTo, val),
                                            trace = FALSE, verbose = FALSE,
                                            start = fittedCoef))
                if (is.null(updated))
                    break
                updatedDev[(sgn + 1)/2 + 1] <- deviance(updated)
                prof[[par]][maxsteps + 1 + sgn * sub] <-
                    sgn * sqrt((deviance(updated) - fittedDev)/disp)
                par.vals[maxsteps + + 1 + sgn * sub,] <- parameters(updated)
            }
            if (all(updatedDev != 0)) {
                quad <- (sum(updatedDev) - 2 * fittedDev)/(2 * margin^2)
                lin <- (fittedDev - updatedDev[1])/margin +
                    quad * (margin - 2 * fittedCoef[i])
                int <- fittedDev - lin * fittedCoef[i] - quad * fittedCoef[i]^2
                ## adjust so roots approx where deviance gives z = zmax
                int.adj <- int - zmax^2 * disp - fittedDev
                for (sgn in c(-1, 1)) {
                    dir <- (sgn + 1)/2 + 1
                    root <- (-lin + sgn * sqrt(lin^2 - 4 * int.adj * quad))/
                        (2 * quad)
                    firstApprox <- par.vals[maxsteps + 1 + sgn * sub, i]
                    ## if likelihood approx quadratic use default stepsize, else
                    if (sgn * (root - firstApprox) > stepsize[dir]) {
                        ## not gone out far enough, check for asymptote
                        val <- fittedCoef[i] + sgn * 1000
                        updated <-
                            suppressWarnings(update(fitted, constrain =
                                                    c(fittedConstrain, i),
                                                    constrainTo =
                                                    c(fittedConstrainTo, val),
                                                    trace = FALSE,
                                                    verbose = FALSE,
                                                    start = fittedCoef))
                        if (!is.null(updated) &&
                            sqrt((deviance(updated) - fittedDev)/disp) < zmax)
                            asymptote[dir] <- TRUE
                    }
                    if (abs(root - firstApprox) > stepsize[dir] &&
                        !asymptote[dir]) {
                        prof[[par]][maxsteps + 1 + sgn * sub] <- 0
                        par.vals[maxsteps + 1 + sgn * sub, ] <- NA
                        stepsize[dir] <- abs(root - fittedCoef[i])/(maxsteps/2)
                    }
                }
            }
        }
        for (sgn in c(-1, 1)) {
            if (trace)
                prattle("\nParameter:", par, c("down", "up")[(sgn + 1)/2 + 1],
                        "\n")
            step <- 0
            init <- parameters(fitted)
            while ((step <- step + 1) <= maxsteps) {
                if (step > 2 &&
                    abs(prof[[par]][maxsteps + 1 + sgn * (step - 2)]) > zmax)
                    break
                if (prof[[par]][maxsteps + 1 + sgn * step] != 0)
                    next
                val <- fittedCoef[i] + sgn * step * stepsize[(sgn + 1)/2 + 1]
                updated <-
                    suppressWarnings(update(fitted, constrain =
                                            c(fittedConstrain, i),
                                            constrainTo =
                                            c(fittedConstrainTo, val),
                                            trace = FALSE, verbose = FALSE,
                                            start = init))
                if (is.null(updated)) {
                    message("Could not complete profile for", par, "\n")
                    break
                }
                init <- parameters(updated)
                zz <- (deviance(updated) - fittedDev)/disp
                if (zz > -0.001)
                  zz <- max(zz, 0)
                else stop("profiling has found a better solution, ",
                          "so original fit had not converged")
                prof[[par]][maxsteps + 1 + sgn * step] <- sgn * sqrt(zz)
                par.vals[maxsteps + 1 + sgn * step,] <- init
                #print(data.frame(step = step, val = bi, deviance = fm$deviance,
                                 #zstat = z))
            }
        }
        prof[[par]] <- structure(data.frame(prof[[par]][!is.na(par.vals[,1])]),
                                 names = "z")
        prof[[par]]$par.vals <- par.vals[!is.na(par.vals[,1]),]
        attr(prof[[par]], "asymptote") <- asymptote
    }
    val <- structure(prof, original.fit = fitted, summary = summ)
    class(val) <- c("profile.gnm", "profile.glm", "profile")
    val
}
