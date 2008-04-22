boxcox.gnm <- function (object, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("boxcox is not implemented for gnm objects")
    else
        NextMethod
}
