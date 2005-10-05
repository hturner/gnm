termPredictors.gnm <- function(object, ...) {
    if (is.null(object$termPredictors)){
        object <- update(object, formula = formula(object),
                         termPredictors = TRUE, start = coef(object),
                         verbose = FALSE, trace = FALSE)
    }
    object$termPredictors
}
