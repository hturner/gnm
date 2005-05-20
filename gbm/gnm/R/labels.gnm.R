labels.gnm <- function(object, ...) {
    labelList <- as.list(attr(gnmTerms(object, object$call$eliminate),
                              "termLabels"))
    termAssign <- attr(model.matrix(object), "assign")[!object$constrain]
    if (object$eliminate) termAssign <- termAssign[-seq(object$eliminate)]
    unique(unlist(labelList)[termAssign])
}

