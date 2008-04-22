variable.names.gnm <- function(object, full = FALSE, ...) {
    if (full)
        names(coef(object))
    else {
        setToZero <- object$constrain[object$constrainTo == 0]
        names(coef(object)[-setToZero])
    }
}
