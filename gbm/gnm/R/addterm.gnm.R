addterm.gnm <- dropterm.gnm <- function (object, scope, scale = 0, test = c("none", "Chisq", 
    "F"), k = 2, sorted = FALSE, trace = FALSE, ...) {
    if (inherits(object, "gnm", TRUE) == 1)
        stop(gsub(".gnm", "", as.character(match.call())[[1]]),
             " is not implemented for gnm objects")
    else
        NextMethod
}
