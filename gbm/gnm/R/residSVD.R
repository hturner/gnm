residSVD <- function(model, fac1, fac2, d = 1) {
    if (!is.null(model$call$data)) {
        Data <- eval(model$call$data, parent.frame())
        fac1 <- eval(match.call()$fac1, Data)
        fac2 <- eval(match.call()$fac2, Data)
    }
    fdata <- data.frame(fac1, fac2)
    if (!is.null(model$na.action)) fdata <- fdata[-model$na.action,]
    X <- data.frame(rw = model$residuals * model$weights, w = model$weights)
    X <- lapply(X, tapply, fdata, sum, simplify = TRUE)
    X <- X$rw/X$w
    X <- svd(naToZero(X), d, d)
    init <- c(t(sqrt(X$d[seq(d)]) * t(X$u)), t(sqrt(X$d[seq(d)]) * t(X$v)))
    names(init) <- rep(c(paste("fac1", levels(fac1), sep = "."),
                         paste("fac2", levels(fac2), sep = ".")), d)
    init
}
