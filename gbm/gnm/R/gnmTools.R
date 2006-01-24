"gnmTools" <-
    function(gnmTerms, gnmData, x, termPredictors)
{
    labelList <- attr(gnmTerms, "parsedLabels")
    prefixList <- attr(gnmTerms, "prefixLabels")
    offsetList <- lapply(attr(gnmTerms, "offset"), function(x) {
        if (!is.null(x))
            eval(parse(text = x), envir = gnmData)
        else
            0
    })

    nFactor <- length(labelList)
    termTools <- factorAssign <- vector(mode = "list", length = nFactor)
    for (i in seq(labelList)) {
        if (inherits(labelList[[i]], "Nonlin")) {
            termTools[[i]] <- eval(attr(labelList[[i]], "call"), gnmData,
                                   environment(gnmTerms))
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
    nTheta <- length(factorAssign)
    nr <- dim(gnmData)[1]
    tmp <- seq(factorAssign) * nr
    first <- c(0, tmp[-nTheta])
    last <- tmp - 1
    nc <- tabulate(factorAssign)
    tmp <- cumsum(nc)
    a <- c(1, tmp[-nFactor] + 1)
    z <- tmp
    storage.mode(first) <- storage.mode(last) <- storage.mode(a) <-
        storage.mode(z) <- "integer"
    X <- baseMatrix <- matrix(1, nrow = nr, ncol = nTheta)
    for (i in seq(termTools)) 
        if (is.matrix(termTools[[i]]))
            baseMatrix[, factorAssign == i] <- termTools[[i]]
    colnames(X) <- names(factorAssign)
    storage.mode(X) <- storage.mode(baseMatrix) <- "double"

    multIndex <- gsub("\.Factor[0-9]+\.", "", unlist(prefixList))
    multIndex[multIndex == ""] <- seq(sum(multIndex == ""))
    
    classID <- sapply(labelList, class)
    plugInStart <- !sapply(lapply(termTools, function(x) x$start), is.null)
    thetaClassID <- classID
    thetaClassID[plugInStart] <- "plugInStart"
    thetaClassID <- structure(thetaClassID[factorAssign],
                                names = names(factorAssign))

    if (classID[1] == "Linear")
        X[, factorAssign == 1] <- termTools[[1]]
    
    if (classID[1] == "Linear" || x || termPredictors) {
        termAssign <- unclass(as.factor(multIndex))[factorAssign]
        if (!is.null(attr(gnmTerms, "termsID")))
            termAssign <- attr(gnmTerms, "termsID")[termAssign]
        if (classID[1] == "Linear") {
            linearAssign <- attr(termTools[[1]], "assign")
            termAssign <- termAssign + max(linearAssign) - 1
            termAssign[thetaClassID == "Linear"] <- linearAssign
        }
    }

    start <- function(scale = 0.1) {
        theta <- structure(runif(nTheta, -1, 1) * scale,
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
                switch(classID[i],
                       "Exp" = exp(.Call("submatprod", baseMatrix,
                       parameterList[[i]], first[a[i]], nr, nc[i])),
                       "Nonlin" = termTools[[i]]$predictor(parameterList[[i]]),
                       .Call("submatprod", baseMatrix, parameterList[[i]],
                             first[a[i]], nr, nc[i], PACKAGE = "gnm",
                             NAOK = TRUE))
        }
        factorList <- mapply("+", factorList, offsetList, SIMPLIFY = FALSE)
        if (term && classID[[1]] == "Linear")
            factorList[[1]] <- t(rowsum(t(termTools[[1]]) * parameterList[[1]],
                                        linearAssign))
        unlistOneLevel(factorList)
    }
    
    predictor <- function(factorList, term = FALSE) {
        termPredictors <- lapply(split(factorList, multIndex), do.call,
                                 what = pprod)
        if (term) {
            if (!is.null(attr(gnmTerms, "termsID")))
                 termPredictors <- lapply(split(termPredictors,
                                          attr(gnmTerms, "termsID")), do.call,
                                          what = psum)
            termPredictors <- do.call("cbind", termPredictors)
            colnames(termPredictors) <-
                c("(Intercept)"[attr(attr(gnmTerms, "terms"), "intercept")],
                  attr(attr(gnmTerms, "terms"), "term.labels"))
        }
        else termPredictors <- rowSums(do.call("cbind", termPredictors))
        termPredictors
    }
    
    localDesignFunction <- function(theta, factorList, ind = NULL) {
        if (!is.null(ind)) {
            a <- ind
            if (factorAssign[ind] > 1)
                ind <- ind - z[factorAssign[ind] - 1]
        }
            
        for (i1 in a) {
            fi <- factorAssign[i1]
            i2 <- ifelse(is.null(ind), z[fi], i1)            
            switch(classID[fi],
                   "Exp" = {
                       v <- do.call("pprod",
                                    factorList[multIndex == multIndex[fi]])
                   },
                   "Nonlin" = {
                       .Call("nonlin", X, first[i1], last[i2],
                             quote(termTools[[fi]]$localDesignFunction(
                             coef = theta[factorAssign == fi],
                             predictor = factorList[[fi]], ind = ind)),
                             environment(), PACKAGE = "gnm")
                   },
                   "character" = {
                       v <- do.call("pprod",
                                    factorList[(multIndex == multIndex[fi])
                                               & (seq(multIndex) != fi)])
                   },
                   v <- numeric(0))
            if (length(v)) {
                    .Call("subprod", X, baseMatrix, as.double(v),
                          first[i1], last[i2], nr, PACKAGE = "gnm")
            }
        }
        if(!is.null(ind)) X[, a, drop = FALSE]
        else X
    }

    toolList <- list(classID = thetaClassID, start = start,
                     factorList = factorList, predictor = predictor,
                     localDesignFunction = localDesignFunction)
    if (x) toolList$termAssign <- termAssign
    toolList
}
