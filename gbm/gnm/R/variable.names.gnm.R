variable.names.gnm <- function(object, full = FALSE, ...) {
    if (full)
        names(coef(object))
    else if (object$eliminate) {
        if (object$eliminate < length(coef(object)))
            names(coef(object)[!is.na(coef(object)) &
                               seq(coef(object)) > object$eliminate])
        else
            cat("No non-eliminated coefficients\n")
    }
    else
        names(coef(object)[!is.na(coef(object))])
}
