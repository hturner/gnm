model.matrix.gnm <- function(object, ...) {
    if (!"x" %in% names(object)){
        object <- update(object, formula = formula(object), x = TRUE,
                         start = coef(object), verbose = FALSE, trace = FALSE)
    }
    return(object[[match("x", names(object))]])
}
