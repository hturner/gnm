add1.gnm <- function(object, scope, scale = 0, test = c("none", "Chisq", 
    "F"), x = NULL, k = 2, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("add1 is not implemented for gnm objects")
    else
        NextMethod
}
