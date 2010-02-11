termPredictors.gnm <- function(object, ...) {
    if (is.null(object$termPredictors)){
        modelData <- model.frame(object)
        modelTools <- gnmTools(terms(object), modelData)
        nelim <- nlevels(object$eliminate)
        if (nelim) theta <- parameters(object)[-seq_len(nelim)]
        varPredictors <- modelTools$varPredictors(theta)
        termPredictors <- modelTools$predictor(varPredictors, term = TRUE)
        rownames(termPredictors) <- rownames(modelData)
        if (nelim) termPredictors <-
            cbind("(eliminate)" = unname(parameters(object)[object$eliminate]),
                  termPredictors)
        termPredictors
    }
    else object$termPredictors
}
