fitted.gnm <- function(object, ...) {
    fitted <- NextMethod("fitted")

    if (!is.null(object$table.attr))
        attributes(fitted) <- object$table.attr

    fitted
}
