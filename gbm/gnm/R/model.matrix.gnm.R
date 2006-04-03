model.matrix.gnm <- function(object, coef = NULL, ...) {
    if (!"x" %in% names(object) || !is.null(coef)) {
        xcall <- object$call
        xcall$method <- "model.matrix"
        xcall$constrain <- object$constrain
        xcall$constrainTo <- object$constrainTo
        xcall$verbose <- xcall$trace <- FALSE
        if (!is.null(coef))
            xcall$start <- coef
        else
            xcall$start <- coef(object)
        extras <- match.call(gnm, expand.dots = FALSE)$...
        if (length(extras) > 0) {
            existing <- !is.na(match(names(extras), names(xcall)))
            for (a in names(extras)[existing]) xcall[[a]] <- extras[[a]]
            if (any(!existing)) {
                xcall <- c(as.list(xcall), extras[!existing])
                xcall <- as.call(xcall)
            }
        }
        env <- environment(formula(object))
        if (is.null(env)) 
            env <- parent.frame()
        Data <- eval(xcall$data, env)   
        eval(xcall, as.data.frame(Data), env)
    }
    else object[[match("x", names(object))]]
}
