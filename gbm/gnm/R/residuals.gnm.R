residuals.gnm <- function(object, type = "deviance", ...) {
    if (type == "partial")
        stop("type = \"partial\" not implemented for gnm objects.")
    else
        NextMethod("residuals")
}
