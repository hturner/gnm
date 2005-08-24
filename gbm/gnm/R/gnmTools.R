"gnmTools" <-
    function(gnmEnvironment, gnmTerms, gnmData, x, termPredictors)
{
    labelList <- attr(gnmTerms, "parsedLabels")
    prefixList <- attr(gnmTerms, "prefixLabels")
    offsetList <- lapply(attr(gnmTerms, "offset"), function(x) {
        if (!is.null(x))
            eval(parse(text = x), envir = gnmData)
        else
            0
    })

    termTools <- factorAssign <- labelList
    for (i in seq(labelList)) {
        if (inherits(labelList[[i]], "Nonlin")) {
            termTools[[i]] <- eval(attr(labelList[[i]], "call"),
                                   envir = gnmData, enclos = gnmEnvironment)
            factorAssign[[i]] <-
                structure(rep(i, length(termTools[[i]]$labels)),
                          names = paste(prefixList[[i]], ".",
                          termTools[[i]]$labels, sep = ""))
        }
        else {
            termTools[[i]] <- model.matrix(terms(reformulate(labelList[[i]]),
                                                 keep.order = TRUE),
                                           data = gnmData)
            factorAssign[[i]] <- structure(rep(i, ncol(termTools[[i]])),
                                           names = paste(prefixList[[i]],
                                           colnames(termTools[[i]]), sep = ""))
        }
    }

    factorAssign <- do.call("c", factorAssign)

    multIndex <- gsub("\.Factor[0-9]+\.", "", unlist(prefixList))
    multIndex[multIndex == ""] <- seq(sum(multIndex == ""))
    
    classID <- sapply(labelList, class)
    plugInStart <- !sapply(lapply(termTools, function(x) x$start), is.null)
    thetaClassID <- classID
    thetaClassID[plugInStart] <- "plugInStart"
    thetaClassID <- structure(thetaClassID[factorAssign],
                                names = names(factorAssign))
    
    if (x | termPredictors) {
        termAssign <- unclass(as.factor(multIndex))[factorAssign]
        if ("Linear" %in% classID) {
            linearAssign <- attr(termTools[[1]], "assign")
            termAssign <- termAssign + max(linearAssign) - 1
            termAssign[thetaClassID == "Linear"] <- linearAssign
        }
    }

    start <- function(scale = 0.1) {
        theta <- structure(runif(length(factorAssign), -1, 1) * scale,
                           names = names(factorAssign))
        theta <- ifelse(theta < 0, theta - scale, theta + scale)
        for (i in seq(termTools)[plugInStart])
            theta[factorAssign == i] <- termTools[[i]]$start
        theta
    }

    factorList <- function(theta, term = FALSE) {
        factorList <- parameterList <- unname(split(theta, factorAssign))
        for (i in seq(factorList)) {
            factorList[[i]] <-
                switch(classID[[i]],
                       "Exp" = exp(drop(termTools[[i]] %*%
                       parameterList[[i]])),
                       "Nonlin" = termTools[[i]]$predictor(parameterList[[i]]),
                       drop(termTools[[i]] %*% parameterList[[i]]))
        }
        factorList <- mapply("+", factorList, offsetList, SIMPLIFY = FALSE)
        if (term & classID[[1]] == "Linear")
            factorList[[1]] <- t(rowsum(t(termTools[[1]] %*%
                                          diag(parameterList[[1]])),
                                        linearAssign))
        unlistOneLevel(factorList)
    }
    
    predictor <- function(factorList, term = FALSE) {
        termPredictors <-
            do.call("cbind", tapply(structure(factorList, class = "list"),
                                    multIndex,
                                    function(list) do.call("pprod", list)))
        if (term) colnames(termPredictors) <-
            c("(Intercept)"[attr(attr(gnmTerms, "terms"), "intercept")],
              attr(gnmTerms, "termLabels"))
        else termPredictors <- rowSums(termPredictors)
        termPredictors
    }
    
    localDesignFunction <- function(theta, factorList) {
        derivativeList <- productList <- list()
        for (i in seq(termTools)) derivativeList[[i]] <- 
            switch(classID[[i]],
                   "Exp" = factorList[[i]] * termTools[[i]],
                   "Nonlin" = termTools[[i]]$localDesignFunction(
                   coef = theta[factorAssign == i],
                   predictor = factorList[[i]]),
                   termTools[[i]])
        for (i in seq(derivativeList)) 
            productList[[i]] <- derivativeList[[i]] * 
                do.call("pprod", (factorList[(multIndex == multIndex[i]) &
                                   (seq(multIndex) != i)]))
        structure(do.call("cbind", productList),
                  dimnames = list(NULL, names(factorAssign)))
    }

    theta <- start()
    unestimable <- apply(localDesignFunction(theta, factorList(theta)) == 0,
                                             2, all)

    toolList <- list(unestimable = unestimable, classID = thetaClassID,
                     start = start, factorList = factorList,
                     predictor = predictor,
                     localDesignFunction = localDesignFunction)
    if (x) toolList$termAssign <- termAssign
    toolList
}
