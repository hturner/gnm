labels.gnm <- function(object, ...) {
    labels <- attr(gnmTerms(formula(object), object$call$eliminate),
                   "termLabels")
    termAssign <- attr(model.matrix(object), "assign")[!object$constrain]
    if (object$eliminate) termAssign <- termAssign[-seq(object$eliminate)]
    unique(labels[termAssign])
}

