kappa.gnm <- function (z, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("kappa is not implemented for gnm objects")
    else
        NextMethod
}    
