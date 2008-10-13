psum <- function(...) {
    summandList <- list(...)
    nSummands <- length(summandList)
    if (nSummands == 0) return(0)
    else if (nSummands == 1) return(summandList[[1]])
    else {
        trySum <- try(summandList[[1]] + do.call("Recall", summandList[-1]),
                      silent = TRUE)
        if (inherits(trySum, "try-error"))
            stop("addition not implemented for types of argument supplied")
        else trySum
    }
}
