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
            termTools[[i]] <- model.matrix(reformulate(labelList[[i]]),
                                           data = gnmData)
            factorAssign[[i]] <- structure(rep(i, ncol(termTools[[i]])),
                                           names = paste(prefixList[[i]],
                                           colnames(termTools[[i]]), sep = ""))
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
    thetaClassIndex <- structure(classIndex[factorAssign],
                                names = names(factorAssign))

    start <- function (scale = 0.1) {
        theta <- structure(runif(length(factorAssign), -1, 1) * scale,
                           names = names(factorAssign))
        theta <- ifelse(theta < 0, theta - scale, theta + scale)
        ind <- thetaClassIndex == "Linear"
        if (any(ind))
            theta[ind] <- naToZero(glm.fit(termTools[[1]],
                                           model.response(gnmData),
                                           weights = weights, offset = offset,
                                           family = family)$coefficients)
        for (i in seq(termTools)[!sapply(lapply(termTools, function(x) x$start),
                                         is.null)]) {
            ind <- factorAssign == i
            thetaClassIndex[ind] <- "plugInStart"
            if (is.function(termTools[[i]]$start))
                theta[ind] <- termTools[[i]]$start(sum(ind))
            else
                theta[ind] <- termTools[[i]]$start
        }
        theta
    }

    factorList <- function(theta) {
        factorList <- parameterList <- unname(split(theta, factorAssign))
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
    
    localDesignFunction <- function(theta, factorList) {
        derivativeList <- productList <- list()
        for (i in seq(termTools)) derivativeList[[i]] <- 
            switch(classIndex[[i]],
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
                  dimnames = list(NULL, colnames = names(factorAssign)))
    }

    toolList <- list(classIndex = thetaClassIndex, start = start,
                     factorList = factorList, predictor = predictor,
                     localDesignFunction = localDesignFunction)
    if (x) toolList$termAssign <- termAssign
    toolList
}
