variable.names.gnm <- function(object, full = FALSE, ...) {
    if (full)
        names(coef(object))[!object$auxiliary]
    else
        names(coef(object)[!is.na(coef(object)) & !object$auxiliary])
}
