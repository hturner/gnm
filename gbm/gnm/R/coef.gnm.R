coef.gnm <- function(object, ...) {
    if (object$eliminate)
        object$coefficients[-seq(object$eliminate)]
    else   
        object$coefficients
}
