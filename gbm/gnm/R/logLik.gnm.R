logLik.gnm <- function (object, ...) 
{
    if (length(list(...))) 
        warning("extra arguments discarded")
    fam <- family(object)$family
    p <- object$rank
    if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) 
        p <- p + 1
    val <- p - object$aic/2 + object$eliminate
    attr(val, "df") <- p - object$eliminate
    class(val) <- "logLik"
    val
}
