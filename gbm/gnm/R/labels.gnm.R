labels.gnm <- function(object) {
    labelList <- as.list(attr(object$terms, "term.labels"))
    termAssign <- attr(model.matrix(object), "assign")[!object$constrain]
    unique(unlist(labelList)[termAssign])
}

