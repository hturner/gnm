model.matrix.gnm <- function(object, ...) {
    if (!"x" %in% names(object)) {
        xCall <- object$call
        xCall$method <- "model.matrix"
        xCall$constrain <- object$constrain
        xCall$start <- coef(object)
        xCall$verbose <- xCall$trace <- FALSE
        xCall[[1]] <- as.name("gnm")
        env <- environment(formula(object))
        if (is.null(env)) 
            env <- parent.frame()
        eval(xCall, env)
    }
    else object[[match("x", names(object))]]
}
