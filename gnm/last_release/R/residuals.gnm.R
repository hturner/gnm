residuals.gnm <- function(object, type = "deviance", ...) {
    if (type == "partial")
        stop("type = \"partial\" not implemented for gnm objects.")
    else if (type == "deviance") {
        y <- object$y
        mu <- object$fitted.values
        wts <- object$prior.weights
        if (object$df.res > 0) {
            res <- sqrt(pmax((object$family$dev.resids)(y, mu, wts), 0))
            res[y < mu] <- -res[y < mu]
        }
        else res <- rep.int(0, length(mu))

        if (!is.null(object$na.action))
            res <- naresid(object$na.action, res)
    }
    else
        res <- NextMethod("residuals")

    if (!is.null(object$table.attr))
        attributes(res) <- object$table.attr

    res
}
