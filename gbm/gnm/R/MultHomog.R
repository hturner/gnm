MultHomog <- function(...){
    labelList <- as.character((match.call(expand.dots = FALSE))[[2]])
    gnmData <- getModelFrame()

    designList <- lapply(labelList, function(x) {
        model.matrix(reformulate(x), data = gnmData)
    })

    # assume factors have same levels with same names
    labels <- colnames(designList[[1]])
    
    predictor <- function(coef) {
        do.call("pprod", lapply(designList, "%*%", coef))
    }

    localDesignFunction <- function(coef, ...) {
        productList <- designList
        for (i in seq(designList))
            productList[[i]] <- designList[[i]] * 
                drop(do.call("pprod", lapply(designList[-i], "%*%", coef)))
        do.call("psum", productList)
    }

    list(labels = labels, predictor = predictor,
         localDesignFunction = localDesignFunction)
}
