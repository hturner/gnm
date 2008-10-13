instances <- function(term, instances = 1){
    term <- match.call()$term
    if (!"inst" %in% names(formals(match.fun(term[[1]]))))
        stop(term[[1]], " has no inst argumnt")
    termList <- vector(mode = "list", length = instances)
    for (i in seq(instances)) {
        termList[[i]] <- term
        termList[[i]]$inst <- i
    }
    paste(unlist(termList), collapse = " + ")
}
