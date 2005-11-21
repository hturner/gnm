model.matrix.gnm <- function(object, ...) {
    if (!"x" %in% names(object))
        return(update(object, formula = formula(object),
                      method = "model.matrix", start = coef(object),
                      verbose = FALSE, trace = FALSE))
    object[[match("x", names(object))]]
}
