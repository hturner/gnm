add1.gnm <- addterm.gnm <- alias.gnm <- anova.gnm <- boxcox.gnm <-
    confint.gnm <- dfbeta.gnm <- dfbetas.gnm <- drop1.gnm <-
    dropterm.gnm <- dummy.coef.gnm <- effects.gnm <- kappa.gnm <-
    influence.gnm <- logtrans.gnm <- predict.gnm <- proj.gnm  <-
    rstudent.gnm <- function(...) {
        stop(gsub(".gnm", "", as.character(match.call())[[1]]),
             " is not implemented for gnm objects")
    }
