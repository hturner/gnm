confint.gnm <- function (object, parm, level = 0.95, trace = FALSE, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("confint is not implemented for gnm objects")
    else
        NextMethod
}
