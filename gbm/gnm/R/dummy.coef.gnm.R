dummy.coef.gnm.R <- function (object, use.na = FALSE, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("dummy.coef is not implemented for gnm objects")
    else
        NextMethod
}    
