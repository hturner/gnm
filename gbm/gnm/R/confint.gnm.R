confint.gnm <- function (object, parm = ofInterest(object), level = 0.95,
                         trace = FALSE, ...) 
{
    pnames <- names(coef(object))
    if (missing(parm)) 
        parm <- seq(along = pnames)
    else if (is.character(parm)) 
        parm <- match(parm, pnames, nomatch = 0)
    cat("Waiting for profiling to be done...\n")
    flush.console()
    object <- profile(object, which = parm, alpha = 1 - level, trace = trace)
    confint(object, level = level, ...)
}
