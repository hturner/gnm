coef.gnm <- function(object) {
  structure(object$coefficients, class = "coef.gnm",
            auxiliary = object$auxiliary)
}
