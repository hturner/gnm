term.predictors.default <- function(object, ...) {
    if (is.null(object$term.predictors)){
        X <- model.matrix(object)
        term.predictors <- t(rowsum(t(X %*% diag(coef(object))),
                                    attr(X, "assign")))
        colnames(term.predictors) <- c("(Intercept)"[attr(X, "assign") == 0],
                                       attr(object$terms, "term.labels"))
        term.predictors
    }
    else object$term.predictors
}
