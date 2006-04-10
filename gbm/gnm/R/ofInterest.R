ofInterest <- function(object, ...) {
    if (!is.null(object$ofInterest)) {
        if (!is.null(coef(object))) {
            if (is.matrix(coef(object)))
                rownames(coef(object))[object$ofInterest]
            else
                names(coef(object))[object$ofInterest]
        }
        else
            object$ofInterest
    }
    else
        NULL
}
