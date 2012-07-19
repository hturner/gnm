hatvalues.gnm <- function(model, ...) {
    X <- as(model.matrix(model), "sparseMatrix")
    var <- unclass(vcov(model, with.eliminate = TRUE))
    eliminate <- model$eliminate
    scale <- model$weights/attr(var, "dispersion")
    hat <- rowSums((X %*% var) * X) * scale
    if (!is.null(eliminate)) {
        ## no covElim!
        if (length(model$constrain))
            X <- X[, -model$constrain, drop = FALSE]
        hat <- hat + (2 * rowSums(X * attr(var, "covElim")[eliminate,]) +
                      attr(var, "varElim")[eliminate]) * scale
    }
    hat <- naresid(model$na.action, hat)
    hat[is.na(hat)] <- 0
    hat[hat > 1 - 100 * .Machine$double.eps] <- 1
    if (!is.null(model$table.attr))
        attributes(hat) <- model$table.attr
    hat
}

