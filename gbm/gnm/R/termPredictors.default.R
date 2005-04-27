termPredictors.default <- function(object, ...) {
    if (is.null(object$termPredictors)){
        X <- model.matrix(object)
        termPredictors <- t(rowsum(t(X %*% diag(coef(object))),
                                    attr(X, "assign")))
        colnames(termPredictors) <- c("(Intercept)"[attr(X, "assign") == 0],
                                       attr(object$terms, "term.labels"))
        termPredictors
    }
    else object$termPredictors
}
