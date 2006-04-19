pick <- function(vec, setlabels = NULL,  title = NULL,
                 items.label = "Pick from:", warn = "Nothing chosen",
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
