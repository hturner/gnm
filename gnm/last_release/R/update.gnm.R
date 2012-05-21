update.gnm <- function (object, formula., ..., evaluate = TRUE)
{
    call <- object$call
    if (is.null(call))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) {
        ## update.formula reorders nonlin terms as lin (main effects)
        ## therefore use substitute to keep order
        rhs <- formula.[[length(formula.)]]
        rhs <- do.call(substitute,
                       list(rhs, env = list("." = object$formula[[3]])))
        if (length(formula.) == 3) {
            lhs <- formula.[[2]]
            lhs <- do.call(substitute,
                           list(lhs, env = list("." = object$formula[[2]])))
            call$formula <- call("~", lhs, rhs)
        } else call$formula <- call("~", rhs)
    }
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate)
        eval(call, parent.frame())
    else call
}

