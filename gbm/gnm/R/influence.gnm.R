influence.gnm <- function (model, do.coef = TRUE, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("influence is not implemented for gnm objects")
    else
        NextMethod
}     
