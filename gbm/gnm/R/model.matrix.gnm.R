model.matrix.gnm <- function(object, ...) {
    if (!"x" %in% names(object)) {
        xcall <- object$call
        xcall$method <- "model.matrix"
        xcall$constrain <- object$constrain
        xcall$start <- coef(object)
        xcall$verbose <- xcall$trace <- FALSE
        env <- environment(formula(object))
        if (is.null(env)) 
            env <- parent.frame()
        if (!is.null(xcall$data))
            Data <- eval(xcall$data, env)
        eval(xcall, as.data.frame(Data), env)
    }
    else object[[match("x", names(object))]]
}
