"ofInterest<-" <- function(object, value = NULL) {
    coefNames <- names(coef(object))
    if (identical(value, "pick")) {
        value <-
            relimp:::pickFrom(coefNames,
                              setlabels = "Coefficients of interest",
                              title = "Select coefficients of interest",
                              items.label = "Model coefficients:",
                              edit.setlabels = FALSE)
        value <- unname(unlist(value))
        if(!length(nchar(value))) {
            warning("No subset of coefficients selected ",
                    "- assuming all are of interest. ")
            value <- NULL
        }
        value <- match(value, coefNames)
    }
    if (is.character(value))
        value <- grep(value, coefNames)
    if (!is.null(value) && !any(value %in% seq(coefNames))) 
        stop("One or more replacement value is invalid.")
    object$ofInterest <- value
    object
}
