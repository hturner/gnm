labels.gnm <- function(object) {
    labelList <- as.list(attr(object$terms, "term.labels"))
    termAssign <- attr(model.matrix(object),
                       "assign")[!object$constrain & !object$auxiliary]
    unique(unlist(labelList)[termAssign])
}

