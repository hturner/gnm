hatvalues.gnm <- function(object) {
    X <- model.matrix(object)
    browser()
    hat <- diag(X %*% summary(object)$cov.unscaled %*% t(X) %*%
              diag(object$weights))
    hat <- naresid(object$na.action, hat)
    hat[is.na(hat)] <- 0
    hat
}
   
