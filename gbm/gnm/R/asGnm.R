asGnm <- function(object, ...){
    if (is.null(object))
        return(NULL)
    UseMethod("as.gnm")
}
