drop1.gnm.R <- function (object, scope, scale = 0, test = c("none", "Chisq", 
    "F"), k = 2, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("drop1 is not implemented for gnm objects")
    else
        NextMethod
}    
