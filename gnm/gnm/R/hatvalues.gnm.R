hatvalues.gnm <- function(model, ...) {
    X <- model.matrix(model)
    var <- vcov(model)
    hat <- diag(X %*% (var/attr(var, "dispersion")) %*% t(X) %*%
              diag(c(model$weights)))
    hat <- naresid(model$na.action, hat)
    hat[is.na(hat)] <- 0
    hat[hat > 1 - 100 * .Machine$double.eps] <- 1
    hat
}
   
