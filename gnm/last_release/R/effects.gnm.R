effects.gnm <- function (object, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("effects is not implemented for gnm objects")
    else
        NextMethod
}
