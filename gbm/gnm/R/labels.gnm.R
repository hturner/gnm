labels.gnm <- function(object, ...) {
    labels <- attr(terms(object), "term.labels")
    termAssign <- attr(model.matrix(object), "assign")[!object$constrain]
    if (object$eliminate) termAssign <- termAssign[-seq(object$eliminate)]
    unique(labels[termAssign])
}

