rstudent.gnm <- function (model, infl = influence(model, do.coef = FALSE),
                          ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("rstudent is not implemented for gnm objects")
    else
        NextMethod
}     
