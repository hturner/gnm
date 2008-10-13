weights.gnm <- function(object, type = c("prior", "working"), ...) {
    weights <- NextMethod("weights")

    if (!is.null(object$table.attr))
        attributes(weights) <- object$table.attr

    weights
}
