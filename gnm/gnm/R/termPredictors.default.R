termPredictors.default <- function(object, ...) {
    if (is.null(object$termPredictors)){
        X <- model.matrix(object)
        termPredictors <- t(rowsum(t(X %*% diag(naToZero(coef(object)))),
                                    attr(X, "assign")))
        colnames(termPredictors) <- c("(Intercept)"[0 %in% attr(X, "assign")],
                                       attr(object$terms, "term.labels"))
        termPredictors
    }
    else object$termPredictors
}
