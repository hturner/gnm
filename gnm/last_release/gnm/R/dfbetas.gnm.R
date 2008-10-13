dfbetas.gnm <- function (model, ...) {
    if (inherits(model, "gnm", TRUE) == 1)
        stop("dfbetas is not implemented for gnm objects")
    else
        NextMethod
}
