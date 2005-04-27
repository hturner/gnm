termPredictors.gnm <- function(object, ...) {
    if (is.null(object$termPredictors)){
        start <- coef(object)
        object <- update(object, termPredictors = TRUE, start = start,
                         trace = FALSE)
    }
    object$termPredictors
}
