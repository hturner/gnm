"gnmTools" <-
    function(modelTerms, gnmData, x, termPredictors)
{
    unitLabels <- attr(modelTerms, "unitLabels")
    common <- attr(modelTerms, "common")
    prefixLabels <- attr(modelTerms, "prefixLabels")
    match <- attr(modelTerms, "match")
    varLabels <- attr(modelTerms, "varLabels")
    block <- attr(modelTerms, "block")
    classID <- attr(modelTerms, "classID")
    NonlinID <- attr(modelTerms, "NonlinID")

    nFactor <- length(varLabels)
    seqFactor <- seq(nFactor)
    termTools <- factorAssign <- thetaID <- vector(mode = "list",
                                                   length = nFactor)
    blockID <- unique(block)
    adj <- 1
    for (i in blockID) {
        b <- block == i
        if (sum(b) == 1 && length(grep("Nonlin\\(", unitLabels[b]))) {
            i <- which(b)
            termTools[[i]] <- eval(parse(text = unitLabels[b])[[1]][[2]],
                                   gnmData, environment(modelTerms))
            nTheta <- length(termTools[[i]]$labels)
            factorAssign[[i]] <-
                structure(rep.int(i, nTheta),
                          names = paste(prefixLabels[i], termTools[[i]]$labels,
                                        sep = ""))
            adj <- adj + nTheta
        }
        else if (all(common[b])) {
            designList <- lapply(unitLabels[b],
                                 function(x) class.ind(gnmData[[x]]))

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
            i <- which(b)
            nm <- paste(prefixLabels[i], labels, sep = "")
            termTools[b] <- designList
            factorAssign[b] <- lapply(i, function(x, nLevels, nm)
                                      structure(rep(x, nLevels), names = nm),
                                      nLevels, nm)
            adj <- adj + nLevels
        }
        else {
            tmp <- model.matrix(terms(reformulate(c(0, unitLabels[b])),
                                      keep.order = TRUE), data = gnmData)
            tmpAssign <- attr(tmp, "assign") + !attr(tmp, "assign")[1]
            tmpAssign <- which(b)[tmpAssign]
            nm <- paste(prefixLabels[tmpAssign], colnames(tmp)[match[b] | sum(b) > 1],
                        sep = "")
            names(tmpAssign) <- nm
            termTools[b] <- lapply(split(1:ncol(tmp), tmpAssign),
                                   function(i, M) M[, i , drop = FALSE], tmp)
            factorAssign[b] <- split(tmpAssign, tmpAssign)
            adj <- adj + length(tmpAssign)
        }
    }

    factorAssign <- unlist(factorAssign)
    uniq <- !(duplicated(block) & common)[factorAssign]
    parLabels <- names(factorAssign)
    nTheta <- length(factorAssign)
    thetaID <- numeric(nTheta)
    thetaID[uniq] <- seq(sum(uniq))
    thetaID[!uniq] <- thetaID[common[factorAssign] & uniq]
    nr <- dim(gnmData)[1]
    tmp <- seq(factorAssign) * nr
    first <- c(0, tmp[-nTheta])
    firstX <- first[thetaID]
    last <- tmp - 1
    lastX <- last[thetaID] + 1
    nc <- tabulate(factorAssign)
    tmp <- cumsum(nc)
    a <- c(1, tmp[-nFactor] + 1)
    z <- tmp
    lt <- last[z] - first[a] + 1
    storage.mode(firstX) <- storage.mode(first) <- storage.mode(lastX) <-
        storage.mode(last) <- storage.mode(a) <-
        storage.mode(z) <- storage.mode(lt) <- "integer"
    baseMatrix <- matrix(1, nrow = nr, ncol = nTheta)
    for (i in seq(termTools))
        if (is.matrix(termTools[[i]]))
            baseMatrix[, factorAssign == i] <- termTools[[i]]
    X <- baseMatrix
    colID <- match(thetaID, thetaID)
    thetaID <- split(thetaID, factorAssign)
    names(thetaID) <- varLabels
    if (any(duplicated(parLabels[uniq]))){
        parLabels[uniq] <- make.unique(parLabels[uniq])
        warning("Using make.unique() to make default parameter labels unique",
                call. = FALSE)
    }
    colnames(X) <- parLabels
    X <- X[, uniq, drop = FALSE]

    thetaClassID <- structure(classID[factorAssign], names = parLabels)
    theta <- rep(NA, nTheta)
    for (i in blockID) {
        b <- block == i
        if (sum(b) == 1 && is.list(termTools[[which(b)]]) &&
            !is.null(termTools[[which(b)]]$start)){
            theta[unlist(thetaID[b])] <- termTools[[which(b)]]$start
            thetaClassID[unlist(thetaID[b])][is.na(termTools[[which(b)]]$start)] <-
                "Linear"
        }
    }
    names(theta) <- parLabels

    thetaClassID <- thetaClassID[uniq]
    theta <- theta[uniq]

    for (i in seq(attr(modelTerms, "predictor"))) {
        if (!is.null(attr(modelTerms, "start")[[i]])) {
            termID <- unique(unlist(thetaID[attr(modelTerms, "assign") == i]))
            theta[termID] <- attr(modelTerms, "start")[[i]](theta[termID])
        }
    }

    if (x || termPredictors) {
        termAssign <- attr(modelTerms, "assign")[factorAssign]
        if (attr(modelTerms, "intercept"))
            termAssign <- termAssign - 1
    }

    prodList <- vector(mode = "list", length = nFactor)
    names(prodList) <- varLabels
    NonlinID <- NonlinID == "Nonlin"
    classID <- classID == "Special"

    varPredictors <- function(theta) {
        for (i in seqFactor) {
            if (NonlinID[i])
                prodList[[i]] <-
                    drop(termTools[[i]]$predictor(theta[thetaID[[i]]]))
            else
                prodList[[i]] <- .Call("submatprod", baseMatrix,
                                       theta[thetaID[[i]]],
                                       first[a[i]], nr, nc[i],
                                       PACKAGE = "gnm", NAOK = TRUE)
        }
        prodList
    }

    predictor <- function(varPredictors, term = FALSE) {
        if (term)
            sapply(attr(modelTerms, "predictor"), eval,
                   varPredictors)
        else
            eval(e, varPredictors)
    }

    gnmData <- lapply(gnmData[, !names(gnmData) %in% varLabels, drop = FALSE],
                      drop)
    e <- do.call("bquote",
                 list(sumExpression(attr(modelTerms, "predictor")), gnmData))
    varDerivs <- lapply(varLabels, deriv, expr = e)


    commonAssign <- factorAssign[colID]
    nCommon <- table(commonAssign[!duplicated(factorAssign)])
    tmpID <- unique(commonAssign)
    tmpID <- tmpID[(NonlinID | classID)[tmpID]]
    nCommon <- as.integer(nCommon[as.character(tmpID)])
    if (any(NonlinID | classID))
        specialVarDerivs <- deriv(e, varLabels[(NonlinID | classID)])
    convID <- colID[uniq]
    vID <- cumsum(c(1, nCommon))[seq(nCommon)]

    localDesignFunction <- function(theta, varPredictors, ind = NULL) {
        if (!any(common)) {
            if (!is.null(ind)){
                i1 <- convID[ind]
                tmpID <- commonAssign[i1]
            }

            for (i in tmpID) {
                fi <- unique(factorAssign[commonAssign == i])
                if (is.null(ind)){
                    i1  <- a[fi][1]
                    i2 <- z[fi][1]
                }
                else {
                    i2 <- i1
                    a <- ind
                    if (factorAssign[ind] > 1)
                        ind <- ind - z[factorAssign[ind] - 1]
                }
                if (NonlinID[fi]) {
                    .Call("nonlin", X, first[i1], last[i2],
                          quote(termTools[[fi]]$localDesignFunction(
                          coef = theta[factorAssign == fi],
                          predictor = varPredictors[[fi]], ind = ind)),
                          environment(), PACKAGE = "gnm")
                    if (classID[fi]) {
                        v <- attr(eval(varDerivs[[fi]], varPredictors),
                                  "gradient")
                        .Call("subprod", X, X, as.double(v),
                              first[i1], last[i2], nr, PACKAGE = "gnm")
                    }
                }
                else if (classID[fi]) {
                    v <- attr(eval(varDerivs[[fi]], varPredictors),
                              "gradient")
                    .Call("subprod", X, baseMatrix, as.double(v),
                          first[i1], last[i2], nr, PACKAGE = "gnm")
                }
            }
            if(!is.null(ind)) X[, a, drop = FALSE]
            else X
        }
        else {
            if (is.null(ind)){
                v <- attr(eval(specialVarDerivs, varPredictors),
                          "gradient")
                .Call("newsubprod", baseMatrix, as.double(v), X,
                      first[a[tmpID]], first[vID], firstX[a[tmpID]],
                      as.integer(length(nCommon)), lt[tmpID], lastX[z[tmpID]],
                      nr, nCommon, max(nCommon), PACKAGE = "gnm")
            }
            else {
                i1 <- convID[ind]
                fi <- unique(factorAssign[commonAssign == commonAssign[i1]])
                v <- list()
                for(j in fi)
                    v[[j]] <- attr(eval(varDerivs[[j]], varPredictors),
                                   "gradient")
                .Call("single", baseMatrix, as.double(unlist(v[fi])),
                      first[i1], lt[fi[1]], nr, as.integer(length(fi)),
                      PACKAGE = "gnm")
            }
        }
    }

    toolList <- list(classID = thetaClassID, start = theta,
                     varPredictors = varPredictors, predictor = predictor,
                     localDesignFunction = localDesignFunction)
    if (x) toolList$termAssign <- termAssign
    toolList
}
