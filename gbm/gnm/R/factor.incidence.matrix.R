factor.incidence.matrix <- function(thefactor) {
    lev <- levels(thefactor)
    nlevels <- length(lev)
    n <- length(thefactor)
    result <- matrix(0, n, nlevels)
    colnames(result) <- lev
    rownames(result) <- names(thefactor)
    for (i in 1:n) result[i, thefactor[i]] <- 1
    result
#    model.matrix(~ -1 + thefactor)
}
    
