anova.gnm <- function (object, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("anova is not implemented for gnm objects")
    else
        NextMethod
}
