hatvalues.gnm <- function(model, ...) {
    X <- model.matrix(model)
    hat <- diag(X %*% summary(model)$cov.unscaled %*% t(X) %*%
              diag(as.vector(model$weights)))
    hat <- naresid(model$na.action, hat)
    hat[is.na(hat)] <- 0
    hat[hat > 1 - 100 * .Machine$double.eps] <- 1
    hat
}
   
