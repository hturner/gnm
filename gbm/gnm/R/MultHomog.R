MultHomog <- function(...){
    labelList <- as.character((match.call(expand.dots = FALSE))[[2]])
    designList <- lapply(list(...), class.ind)

    ## get labels for all levels
    allLevels <- lapply(designList, colnames)
    labels <- unique(unlist(allLevels))
    nLevels <- length(labels)

    ## expand design matrices if necessary
    if (!all(mapply(identical, allLevels, list(labels)))) {
        labels <- sort(labels)
        M <- matrix(0, nrow = nrow(designList[[1]]), ncol = nLevels,
                    dimnames = list(NULL, labels))
        designList <- lapply(designList, function(design, M) {
            M[,colnames(design)] <- design
            M}, M)
    }
    
    predictor <- function(coef) {
        do.call("pprod", lapply(designList, "%*%", coef))
    }

    localDesignFunction <- function(coef, ind = NULL, ...) {
        productList <- designList
        for (i in seq(designList)) {
            if (is.null(ind)) 
                productList[[i]] <- designList[[i]] * 
                    drop(do.call("pprod", lapply(designList[-i], "%*%", coef)))
            else
                productList[[i]] <- designList[[i]][, ind] * 
                    drop(do.call("pprod", lapply(designList[-i], "%*%", coef)))
        }
        do.call("psum", productList)
    }

    list(labels = labels, predictor = predictor,
         localDesignFunction = localDesignFunction)
}
