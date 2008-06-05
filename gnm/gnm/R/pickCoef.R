pickCoef <- function(object, regexpr = NULL, ...){
    object <- names(coef(object))
    if (is.null(object))
        stop("Coefficient names cannot be extracted from 'object'")
    if (missing(regexpr)) {
        default <- list(setlabels = "Selected coefficients",
                        title = "Select coefficients of interest",
                        items.label = "Model coefficients:",
                        return.indices = TRUE, edit.setlabels = FALSE,
                        warningText =  "No subset of coefficients selected")
        dots <- list(...)
        dotArgs <- match(names(default), names(dots))
        allArgs <- c(list(object), dots, default[is.na(dotArgs)])
        selection <- unname(unlist(do.call(relimp::pickFrom, allArgs)))
    }
    else
        selection <- grep(regexpr, object)

    if (!length(selection))
        selection <- NULL
    else
        names(selection) <- object[selection]
    selection
}
