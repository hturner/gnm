pprod <- function(...) {
    factorList <- list(...)
    nFactors <- length(factorList)
    if (nFactors == 0) return(1)
    else if (nFactors == 1) return(factorList[[1]])
    else {
        tryProduct <- try(factorList[[1]] * do.call("Recall", factorList[-1]),
                          silent = TRUE)
        if (inherits(tryProduct, "try-error"))
            stop("multiplication not implemented for types of argument supplied")
        else tryProduct
    }
}

