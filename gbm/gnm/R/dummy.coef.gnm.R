dummy.coef.gnm <- function (object, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("dummy.coef is not implemented for gnm objects")
    else
        NextMethod
}    
