model.matrix.gnm <- function(object, ...) {
    if (!"x" %in% names(object)) {
        xCall <- update(object, formula = formula(object),
                        constrain = object$constrain,
                        method = "model.matrix", start = coef(object),
                        verbose = FALSE, trace = FALSE, evaluate = FALSE)
        env <- environment(formula(object))
        if (is.null(env)) 
            env <- parent.frame()
        eval(xCall, env)
    }
    object[[match("x", names(object))]]
}
