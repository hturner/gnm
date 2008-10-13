alias.gnm <- function (object, ...){
    if (inherits(object, "gnm", TRUE) == 1)
        stop("alias is not implemented for gnm objects")
    else
        NextMethod
}
