labels.gnm <- function(object, ...) {
    labels <- attr(terms(object), "term.labels")
    termAssign <- attr(model.matrix(object), "assign")
    if (length(object$constrain))
        termAssign <- termAssign[-object$constrain]
    unique(labels[termAssign])
}

