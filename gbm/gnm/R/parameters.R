parameters <- function(object)
    replace(coef(object), object$constrain, object$constrainTo)
