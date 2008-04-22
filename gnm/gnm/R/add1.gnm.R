add1.gnm <- function(object, scope, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("add1 is not implemented for gnm objects")
    else
        NextMethod
}
