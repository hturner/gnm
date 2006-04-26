instances <- function(term, instances = 1){
    term <- match.call()$term
    termList <- vector(mode = "list", length = instances)
    for (i in seq(instances)) {
        termList[[i]] <- term
        termList[[i]]$inst <- i
    }
    paste(unlist(termList), collapse = " + ")
}
