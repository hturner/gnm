"update.gnm" <-
    function (object, formula., ..., evaluate = TRUE) 
{
    call <- object$call
    has.original.call <- FALSE
    if (!is.null(object$original.call)) {
        has.original.call <- TRUE
        call <- object$original.call
    }
    extras <- match.call(expand.dots = FALSE)$...  
    if (!missing(formula.)) 
        if (has.original.call) {
            call$formula <- update.formula(call$formula, formula.)
        } else  call$formula <- update.formula(call$formula, formula.)
    if (length(extras) > 0) {
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
