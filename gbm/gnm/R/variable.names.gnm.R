variable.names.gnm <- function(object, full = FALSE, ...) {
    if (full)
        names(coef(object))
    else if (object$eliminate)
        names(coef(object)[!is.na(coef(object)) &
                           seq(coef(object)) > object$eliminate])
    else
        names(coef(object)[!is.na(coef(object))])
}
