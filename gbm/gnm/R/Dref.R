Dref <- function(..., formula = ~ 1) {
    labelList <- as.character((match.call(expand.dots = FALSE))[[2]])
    gnmData <- getModelFrame()
    
    # get design matrices for Dref factors
    designList <- lapply(gnmData[, labelList], class.ind)

    ## get labels for global parameters
    allLevels <- lapply(designList, colnames)
    global <- unique(unlist(allLevels))
    nGlobal <- length(global)

    ## get design matrix for local structure
    localData <- model.frame(formula, data = gnmData)
    local <- model.matrix(formula, data = localData)

    ## create index and labels for parameters
    factorIndex <- c(rep(seq(labelList), each = ncol(local)), rep(0, nGlobal))
    if (ncol(local) > 1)
        labels <- c(as.vector(sapply(labelList, paste, colnames(local),
                                     sep = ".")), global)
    else
        labels <- c(labelList, global)

    ## pad out design matrices for Dref factors if necessary
    if (!all(mapply(identical, allLevels, list(global)))) {
        labels <- sort(labels)
        M <- matrix(0, nrow = nrow(gnmData), ncol = nGlobal,
                    dimnames = list(NULL, global))
        designList <- lapply(designList, function(design, M) {
            M[,colnames(design)] <- design
            M}, M)
    }

    predictor <- function(coef) {
        # calculate constrained weights
        W <- matrix(nrow = nrow(gnmData), ncol = length(designList))
        for (i in seq(ncol(W)))
            W[,i] <- drop(exp(local %*% coef[factorIndex == i]))
        W <- W/rowSums(W)
        gamma <- sapply(designList, "%*%", coef[factorIndex == 0])
        predictor <- W * gamma
        structure(predictor, W = W)
    }

    localDesignFunction <- function(predictor, ...) {
        W <- attr(predictor, "W")
        Dintercept <- predictor - W * rowSums(predictor)
        do.call("cbind", c(tapply(Dintercept, col(Dintercept), "*", local),
                           list(do.call("psum", (mapply("*", split(W, col(W)),
                                                        designList,
                                                        SIMPLIFY = FALSE))))))
    }
    
    list(start = rep(0.5, length(labels)), labels = labels,
         predictor = predictor, localDesignFunction = localDesignFunction)
}
