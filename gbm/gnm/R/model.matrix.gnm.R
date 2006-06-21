model.matrix.gnm <- function(object, coef = NULL, ...) {
    if (!"x" %in% names(object) || !is.null(coef)) {
        xcall <- object$call
        xcall$method <- "model.matrix"
        xcall$constrain <- object$constrain
        xcall$constrainTo <- object$constrainTo
        xcall$data <- model.frame(object)
        xcall[c("weights", "offset")] <- NULL
        xcall$verbose <- FALSE
        if (!is.null(coef))
            xcall$start <- coef
        else
            xcall$start <- coef(object)
        eval(xcall)
    }
    else object[[match("x", names(object))]]
}
