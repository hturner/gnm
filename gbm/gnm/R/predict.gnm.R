predict.gnm <- function (object, newdata = NULL, type = c("link", "response", 
    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL, 
    na.action = na.pass, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("predict is not implemented for gnm objects")
    else
        NextMethod
}     
