dfbeta.gnm <- dfbetas.gnm <- function (model, infl = lm.influence(model,
                                              do.coef = TRUE), ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop(gsub(".gnm", "", as.character(match.call())[[1]]),
             " is not implemented for gnm objects")
    else
        NextMethod
}
