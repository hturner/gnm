Dref <- function(..., formula = ~ 1) {
    labelList <- as.character((match.call(expand.dots = FALSE))[[2]])
    gnmData <- getModelFrame()
    
    # get design matrices for Dref factors
    designList <- lapply(labelList, function(x) {
        M <- model.matrix(reformulate(c(-1, x)), data = gnmData)
        colnames(M) <- with(gnmData, levels(as.factor(eval(parse(text = x)))))
        M
    })
    global <- unique(unlist(lapply(designList, colnames)))
    nGlobal <- length(global)

    ## pad out if necessary
    if (!all(mapply(identical, sapply(designList, colnames), global))) {
        M <- matrix(0, nrow = nrow(rowData), ncol = nGlobal,
                    dimnames = list(NULL, global))
        designList <- lapply(designList, function(design, M) {
            M[,colnames(design)] <- design
            M}, M)
    }

    # get design matrix for local structure
    local <- model.matrix(formula, data = gnmData)

    # create index and labels for parameters
    factorIndex <- c(rep(seq(labelList), ncol(local)), rep(0, nGlobal))
    if (ncol(local) > 1)
        labels <- c(as.vector(sapply(labelList, paste, colnames(local),
                                     sep = ".")), global)
    else
        labels <- c(labelList, global)

    predictor <- function(coef) {
        pList <- list()
        for (i in seq(labelList))
            pList[[i]] <- drop(exp(local %*% coef[factorIndex == i])) *
                designList[[i]]
        predictor <- sapply(pList, "%*%", coef[factorIndex == 0])
        structure(predictor, global = do.call("psum", pList))
    }

    localDesignFunction <- function(predictor, ...) {
        do.call("cbind", c(tapply(predictor, col(predictor), "*", local),
                           list(attr(predictor, "global"))))
    }
    
    list(labels = labels, predictor = predictor,
         localDesignFunction = localDesignFunction)
}
