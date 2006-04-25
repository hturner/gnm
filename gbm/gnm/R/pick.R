pick <- function(object, regexp = NULL, setlabels = "Coefficients of interest",
                     title ="Select coefficients of interest",
                     items.label = "Model coefficients:",
                     warn = "No coefficients were chosen",
                     return.indices = TRUE, edit.setlabels = FALSE, ...){
    selection <- relimp:::pickFrom(vec, setlabels = setlabels, title = title,
                                   items.label = items.label,
                                   return.indices = return.indices,
                                   edit.setlabels = edit.setlabels, ...)
    selection <- unname(unlist(selection))
    if(!length(selection)) {
        warning(warn, call. = FALSE)
        selection <- NULL
    }
    selection
}
