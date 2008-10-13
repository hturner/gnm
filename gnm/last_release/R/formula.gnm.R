formula.gnm <- function(x, ...) {
    form <- formula(x$terms)
    if (x$eliminate > 0) {
        tmp <- .Internal(update.formula(form,
                                        substitute(~ . - e + 1,
                                                   list(e = x$call$eliminate))))
        environment(tmp) <- environment(form)
        form <- formula(terms.formula(tmp, simplify = TRUE, keep.order = TRUE))
    }
    return(form)
}
    
