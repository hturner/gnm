drop1.gnm <- function (object, scope, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("drop1 is not implemented for gnm objects")
    else
        NextMethod
}    
