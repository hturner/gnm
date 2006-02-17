model.matrix.gnm <- function(object, coef = NULL, ...) {
    if (!"x" %in% names(object)) {
        xcall <- object$call
        xcall$method <- "model.matrix"
        xcall$constrain <- object$constrain
        xcall$verbose <- xcall$trace <- FALSE
        if (!is.null(coef))
            xcall$start <- coef
        else
            xcall$start <- coef(object)
        env <- environment(formula(object))
        if (is.null(env)) 
            env <- parent.frame()
        if (!is.null(xcall$data))
            Data <- eval(xcall$data, env)
        eval(xcall, as.data.frame(Data), env)
    }
    else object[[match("x", names(object))]]
}
