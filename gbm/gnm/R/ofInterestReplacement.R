"ofInterest<-" <- function(object, value = NULL) {
    coefNames <- names(coef(object))
    if (identical(value, pick))
        value <- pick(coefNames,
                      setlabels = "Coefficients of interest",
                      title = "Select coefficients of interest",
                      items.label = "Model coefficients:",
                      warn = paste("No subset of coefficients selected",
                     "- assuming all are of interest."))
    if (is.character(value))
        value <- grep(value, coefNames)
    if (!is.null(value)) {
        if (!any(value %in% seq(coefNames))) 
            stop("One or more replacement value is invalid.")
        names(value) <- coefNames[value]
    }
    object$ofInterest <- value
    print(value)
    object
}
