coef.gnm <- function(object, ...) {  
    structure(object$coefficients, eliminate = object$eliminate,
              class = "coef.gnm")
}
