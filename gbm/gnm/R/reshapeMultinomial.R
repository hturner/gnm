reshapeMultinomial <- function(data, catvar, keep.as.ordered = FALSE,
                                countvar = "count", idvar = "id") {
    n <- nrow(data)
    ncat <- nlevels(data[[catvar]])
    tmpID <- gl(n, ncat)
    newData <- data[tmpID, ]
    newData[[catvar]] <- gl(ncat, 1, n * ncat, labels = levels(data[[catvar]]),
                            ordered = keep.as.ordered)
    newData[[idvar]] <- .rowID
    newData[[countvar]] <- as.vector(t(class.ind(data[[catvar]])))
    newData
}
