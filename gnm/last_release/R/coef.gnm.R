coef.gnm <- function(object, ...) {
    structure(object$coefficients, elim.coefs = object$elim.coefs,
              ofInterest = object$ofInterest,
              class = "coef.gnm")
}
