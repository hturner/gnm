coef.gnm <- function(object, ...) {
    structure(object$coefficients, ofInterest = object$ofInterest,
              class = c("coef.gnm", "numeric"))
}
