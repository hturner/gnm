hatvalues.gnm <- function(model, ...) {
    X <- model.matrix(model)
    eliminate <- model$eliminate
    var <- vcov(model)
    scale <- model$weights/attr(var, "dispersion")
    hat <- diag(X %*% tcrossprod(var, X)) * scale
    if (!is.null(eliminate))
        hat <- hat + (2 * rowSums(X * attr(var, "covElim")[eliminate,]) +
        attr(var, "varElim")[eliminate]) * scale
    hat <- naresid(model$na.action, hat)
    hat[is.na(hat)] <- 0
    hat[hat > 1 - 100 * .Machine$double.eps] <- 1
    if (!is.null(model$table.attr))
        attributes(hat) <- model$table.attr
    hat
}

