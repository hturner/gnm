coef.gnm <- function(object, ...) {  
    structure(object$coefficients, ofInterest = object$ofInterest,
              class = "coef.gnm")
}
