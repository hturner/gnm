pickCoef <- function(object, regexpr = NULL, ...){
    nelim <- nlevels(object$eliminate)
    coefs <- coef(object)
    coefs <- names(coefs[(elim + 1):length(coefs)])
    if (is.null(coefs))
        stop("Coefficient names cannot be extracted from 'coefs'")
    if (missing(regexpr)) {
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
    else
        selection <- grep(regexpr, coefs)

    if (!length(selection))
        selection <- NULL
    else
        names(selection) <- coefs[selection]
    selection
}
