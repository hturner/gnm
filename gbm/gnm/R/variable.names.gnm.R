variable.names.gnm <- function(object, full = FALSE, ...) {
    if (full)
        names(coef(object))
    else
        names(coef(object)[!is.na(coef(object))])
}
