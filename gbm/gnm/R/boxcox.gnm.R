boxcox.gnm <- function (object, lambda = seq(-2, 2, 1/10), plotit = TRUE,
                        interp = (plotit && (m < 100)), eps = 1/50,
                        xlab = expression(lambda), ylab = "log-Likelihood", 
                        ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop("boxcox is not implemented for gnm objects")
    else
        NextMethod
}
