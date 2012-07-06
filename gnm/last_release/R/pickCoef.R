pickCoef <- function(object, pattern = NULL, value = FALSE, ...){
    coefs <- names(coef(object))
    if (is.null(coefs))
        stop("Coefficient names cannot be extracted from 'object'")
    if (missing(pattern)) {
        default <- list(setlabels = "Selected coefficients",
                        title = "Select coefficients of interest",
                        items.label = "Model coefficients:",
                        return.indices = TRUE, edit.setlabels = FALSE,
                        warningText =  "No subset of coefficients selected")
        dots <- list(...)
        dotArgs <- match(names(default), names(dots))
        allArgs <- c(list(coefs), dots, default[is.na(dotArgs)])
        selection <- unname(unlist(do.call(relimp::pickFrom, allArgs)))
    }
    else {
        selection <- grep(pattern, coefs, value = FALSE, ...)
    }

    if (!length(selection))
        selection <- NULL
    else if (!value)
        names(selection) <- coefs[selection]
    else
        selection <- parameters(object)[selection]
    selection
}
