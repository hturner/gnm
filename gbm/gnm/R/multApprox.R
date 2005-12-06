multApprox <- function(vec = model$residuals, fac1, fac2, d = 1,
                       weights = model$weights, model = NULL) {
    if (!is.null(model$call$data)) {
        Data <- eval(model$call$data, parent.frame())
        fac1 <- eval(match.call()$fac1, Data)
        fac2 <- eval(match.call()$fac2, Data)
    }
    if (!is.factor(fac1)) stop("fac1 must be a factor")
    if (!is.factor(fac2)) stop("fac2 must be a factor")
    if (!is.null(weights) && (length(vec) != length(weights))) {
        stop("vec and weights have different lengths")
    }
    fdata <- na.omit(data.frame(fac1, fac2))
    browser()
    if (!is.null(model$na.action)) fdata <- fdata[-model$na.action,]
    X <- cbind(vec * weights, weights)
    X <- aggregate(X, fdata, sum, na.rm = TRUE)
    X <- xtabs(V1/V2 ~ fac1 + fac2, X)
    X <- svd(X, d, d)
    result <- c(t(sqrt(X$d[seq(d)]) * t(X$u)), t(sqrt(X$d[seq(d)]) * t(X$v)))
    names(result) <- rep(c(paste("fac1", levels(fac1), sep = "."),
                           paste("fac2", levels(fac2), sep = ".")), d)
    result
}
