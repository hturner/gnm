proj.gnm <- function (object, onedf = FALSE, unweighted.scale = FALSE, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("proj is not implemented for gnm objects")
    else
        NextMethod
}     
