labels.gnm <- function(object, ...) {
    labels <- attr(terms(object), "term.labels")
    termAssign <- attr(model.matrix(object), "assign")
    if (length(object$constrain)) {
        isConstrained <- is.element(names(coef(object)), object$constrain)
        termAssign <- termAssign[!isConstrained]
    }
    unique(labels[termAssign])
}

