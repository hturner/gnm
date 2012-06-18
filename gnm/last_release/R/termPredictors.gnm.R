termPredictors.gnm <- function(object, ...) {
    if (is.null(object$termPredictors)){
        modelData <- model.frame(object)
        modelTerms <- terms(object)
        if (!is.empty.model(modelTerms)) {
            modelTools <- gnmTools(modelTerms, modelData)
            theta <- parameters(object)
            varPredictors <- modelTools$varPredictors(theta)
            termPredictors <- modelTools$predictor(varPredictors, term = TRUE)
            rownames(termPredictors) <- rownames(modelData)
        }
        else termPredictors <- modelData[,0]
        if (!is.null(object$eliminate)) termPredictors <-
            cbind("(eliminate)" =
                  as.vector(attr(coef(object), "eliminated")[object$eliminate]),
                  termPredictors)
        termPredictors
    }
    else object$termPredictors
}
