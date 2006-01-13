reshapeMultinomial <- function(data, catvar, countvar = "count", idvar = "id",
                               keep.as.ordered = FALSE) {
    n <- nrow(data)
    ncat <- nlevels(data[[catvar]])
    caseID <- gl(n, ncat)
    newData <- data[caseID, ]
    newData[[catvar]] <- gl(ncat, 1, n * ncat, labels = levels(data[[catvar]]),
                            ordered = keep.as.ordered)
    newData[[countvar]] <- as.vector(t(class.ind(data[[catvar]])))
    newData[[idvar]] <- caseID
    newData
}
