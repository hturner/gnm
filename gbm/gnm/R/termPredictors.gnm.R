termPredictors.gnm <- function(object, ...) {
    if (is.null(object$termPredictors)){
        modelTools <- gnmTools(terms(object), object$model)
        varPredictors <- modelTools$varPredictors(naToZero(coef(object)))
        modelTools$predictor(varPredictors, term = TRUE)
    }
    else object$termPredictors
}
