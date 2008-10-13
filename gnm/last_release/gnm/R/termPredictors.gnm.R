termPredictors.gnm <- function(object, ...) {
    if (is.null(object$termPredictors)){
        modelData <- model.frame(object)
        modelTools <- gnmTools(terms(object), modelData)
        varPredictors <- modelTools$varPredictors(parameters(object))
        termPredictors <- modelTools$predictor(varPredictors, term = TRUE)
        rownames(termPredictors) <- rownames(modelData)
        termPredictors
    }
    else object$termPredictors
}
