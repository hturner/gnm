dropterm.gnm <- function (object, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("dropterm is not implemented for gnm objects")
    else
        NextMethod
}
