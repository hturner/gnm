influence.gnm <- function (model, do.coef = TRUE, ...) {
    if (inherits(model, "gnm", TRUE) == 1)
        stop("influence is not implemented for gnm objects")
    else
        NextMethod
}     
