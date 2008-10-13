rstudent.gnm <- function (model, ...) {
    if (inherits(model, "gnm", TRUE) == 1)
        stop("rstudent is not implemented for gnm objects")
    else
        NextMethod
}     
