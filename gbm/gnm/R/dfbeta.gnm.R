dfbeta.gnm <- function (model, ...) {
    if (inherits(model, "gnm", TRUE) == 1)
        stop("dfbeta is not implemented for gnm objects")
    else
        NextMethod
}
