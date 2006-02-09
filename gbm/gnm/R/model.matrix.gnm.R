model.matrix.gnm <- function(object, ...) {
    if (!"x" %in% names(object)) {
        args <- object$call
        args$method <- "model.matrix"
        args$constrain <- object$constrain
        args$start <- coef(object)
        args$verbose <- args$trace <- FALSE
        args[[1]] <- as.name("list")
        env <- environment(formula(object))
        if (is.null(env)) 
            env <- parent.frame()
        do.call("gnm", eval(args, env))
    }
    else object[[match("x", names(object))]]
}
