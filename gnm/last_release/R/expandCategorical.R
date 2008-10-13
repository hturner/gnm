expandCategorical <- function(data, catvar, sep = ".", countvar = "count",
                              idvar = "id", as.ordered = FALSE) {
    n <- nrow(data)
    inter <- interaction(data[[catvar]])
    ncat <- nlevels(inter)
    caseID <- gl(n, ncat)
    newData <- data[caseID, -match(catvar, names(data))]
    newData[[levels(interaction(catvar, sep = sep))]] <-
        gl(ncat, 1, n * ncat, labels = levels(inter), ordered = as.ordered)
    newData[[countvar]] <- as.vector(t(class.ind(inter)))
    newData[[idvar]] <- caseID
    newData
}
