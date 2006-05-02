"ofInterest<-" <- function(object, value = NULL) {
    coefNames <- names(coef(object))
    if (!is.null(value)) {
        if (!any(value %in% seq(coefNames))) 
            stop("One or more replacement value is invalid.")
        names(value) <- coefNames[value]
    }
    object$ofInterest <- value
    print(value)
    object
}
