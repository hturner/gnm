labels.gnm <- function(object, ...) {
    labels <- attr(terms(object), "term.labels")
    termAssign <- attr(model.matrix(object), "assign")
    if (!is.null(object$constrain))
        termAssign <- termAssign[-object$constrain[,1]]
    unique(labels[termAssign])
}

