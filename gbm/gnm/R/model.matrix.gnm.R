model.matrix.gnm <- function(object, ...) {
    if (!"x" %in% names(object)){
        start <- coef(object)
        object <- update(object, x = TRUE, start = start, trace = FALSE)
    }
    return(object[[match("x", names(object))]])
}
