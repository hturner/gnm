model.matrix.gnm <- function(object, ...) {
    if (!"x" %in% names(object)){
        object <- update(object, x = TRUE, start = coef(object), trace = FALSE)
    }
    return(object[[match("x", names(object))]])
}
