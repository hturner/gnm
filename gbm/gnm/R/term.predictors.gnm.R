term.predictors.gnm <- function(object, ...) {
    if (is.null(object$term.predictors)){
        start <- coef(object)
        object <- update(object, term.predictors = TRUE, start = start,
                         trace = FALSE)
    }
    object$term.predictors
}
