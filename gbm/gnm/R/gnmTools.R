"gnmTools" <-
    function(gnmTerms, gnmData, x, y, family, weights, offset)
{
    labelList <- attr(gnmTerms, "parsed.labels")
    prefixList <- attr(gnmTerms, "prefix.labels")
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
                                    envir = gnmData)
            factorAssign[[i]] <-
                structure(rep(i, length(termTools[[i]]$labels)),
                          names = paste(prefixList[[i]], ".",
                          termTools[[i]]$labels, sep = ""))
        }
        else {
            termTools[[i]] <-
                model.matrix(as.formula(paste("~ - 1 + ", labelList[[i]],
                                              sep = "")), data = gnmData)
            factorAssign[[i]] <- structure(rep(i, ncol(termTools[[i]])),
                                             names = paste(prefixList[[i]],
                                             colnames(termTools[[i]]),
                                             sep = ""))
        }
    }

    factorAssign <- do.call("c", factorAssign)

    multIndex <- gsub("\.Factor[0-9]+\.", "", unlist(prefixList))
    multIndex[multIndex == ""] <- seq(sum(multIndex == ""))

    if (x) {
        termAssign <- unclass(as.factor(multIndex))
        termAssign <- termAssign[factorAssign] - attr(attr(gnmTerms, "terms"),
                                                     "intercept")
    }
    
    classIndex <- sapply(labelList, class)
    classIndex[is.element(multIndex, multIndex[duplicated(multIndex)]) &
               classIndex != "Exp"] <- "MultNotExp"

    start <- function (scale = 0.2) {
        theta <- structure(runif(length(factorAssign), -1, 1) * scale,
                           names = names(factorAssign))
        for (i in seq(termTools)[!sapply(sapply(termTools, attr, "start"),
                                         is.null)]) {
            ind <- factorAssign == i
            if (is.function(termTools[[i]]$start))
                theta[ind] <- termTools[[i]]$start(sum(ind))
            else
                theta[ind] <- termTools[[i]]$start
        }
        ind <- classIndex[factorAssign] == "MultNotExp"
        theta[ind] <- 2 * scale + theta[ind]
        ind <- classIndex[factorAssign] == "character"
        theta[ind] <- 0
        lin <- model.matrix(reformulate(labelList[classIndex == "character"]),
                            data = gnmData)
        lin <- naToZero(glm.fit(lin, model.response(gnmData), weights = weights,
                                offset = offset, family = family)$coefficients)
        theta[names(lin)] <- lin
        theta
    }

    factorList <- function(parameterList) {
        factorList <- parameterList
        for (i in seq(factorList)) {
            factorList[[i]] <-
                switch(classIndex[[i]],
                       "Exp" = exp(drop(termTools[[i]] %*%
                       parameterList[[i]])),
                       "Nonlin" = termTools[[i]]$predictor(parameterList[[i]]),
                       drop(termTools[[i]] %*% parameterList[[i]]))
        }
        mapply("+", factorList, offsetList, SIMPLIFY = FALSE)
    }
    
    predictor <- function(factorList)
        rowSums(do.call("cbind",
                        tapply(structure(factorList, class = "list"),
                               multIndex,
                               function(list) do.call("pprod", list))))
    
    localDesignFunction <- function(parameterList, factorList) {
        derivativeList <- productList <- list()
        for (i in seq(termTools)) derivativeList[[i]] <- 
            switch(classIndex[[i]],
                   "Exp" = factorList[[i]] * termTools[[i]],
                   "Nonlin" = termTools[[i]]$localDesignFunction(
                   coef = parameterList[[i]], predictor = factorList[[i]]),
                   termTools[[i]])
        for (i in seq(derivativeList)) 
            productList[[i]] <- derivativeList[[i]] * 
                do.call("pprod", (factorList[(multIndex == multIndex[i]) &
                                   (seq(multIndex) != i)]))
        do.call("cbind", productList)
    }

    toolList <- list(factorAssign = factorAssign, start = start,
                     factorList = factorList, predictor = predictor,
                     localDesignFunction = localDesignFunction)
    if (x) toolList$termAssign <- termAssign
    toolList
}
