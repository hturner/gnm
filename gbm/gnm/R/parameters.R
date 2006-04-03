parameters <- function(object)
    replace(coef(object), object$constrain + object$eliminate,
            object$constrainTo)
