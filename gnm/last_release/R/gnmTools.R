"gnmTools" <- function(modelTerms, gnmData = NULL, x = TRUE)
{
    eliminate <- attr(modelTerms, "eliminate")
    unitLabels <- attr(modelTerms, "unitLabels")
    common <- attr(modelTerms, "common")
    prefixLabels <- attr(modelTerms, "prefixLabels")
    match <- attr(modelTerms, "match")
    varLabels <- attr(modelTerms, "varLabels")
    block <- attr(modelTerms, "block")
    type <- attr(modelTerms, "type")

    nFactor <- length(varLabels)
    seqFactor <- seq(nFactor)
    termTools <- factorAssign <- thetaID <- vector(mode = "list",
                                                   length = nFactor)
    blockID <- unique(block)
    adj <- 1
    for (i in blockID) {
        b <- block == i
        if (all(common[b])) {
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
            intercept <- as.numeric(i == 0 && eliminate)
            tmp <- model.matrix(terms(reformulate(c(intercept, unitLabels[b])),
                                      keep.order = TRUE), data = gnmData)
            tmpAssign <- attr(tmp, "assign")
            if (intercept) {
                tmp <- tmp[,-1]
                tmpAssign <- tmpAssign[-1]
            }
            tmpAssign <- which(b)[tmpAssign + !tmpAssign[1]]
            ## paste unless prefix == "" and identical(colnames(tmp), "(Intercept)")
            if (length(tmpAssign) == 1 && tmpAssign == 0)
                nm <- "(Intercept)"
            else
                nm <- paste(prefixLabels[tmpAssign], colnames(tmp), sep = "")
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

    theta <- rep(NA, nTheta)
    for (i in blockID) {
        b <- block == i
        if (sum(b) == 1 && is.list(termTools[[which(b)]]) &&
            !is.null(termTools[[which(b)]]$start)){
            theta[unlist(thetaID[b])] <- termTools[[which(b)]]$start
        }
    }
    names(theta) <- parLabels

    termAssign <- attr(modelTerms, "assign")[factorAssign]
    block <- block[factorAssign]
    for (i in seq(attr(modelTerms, "predictor"))) {
        if (!is.null(attr(modelTerms, "start")[[i]])) {
            termID <- termAssign == i & uniq
            split <- block[termID]
            split <- match(split, unique(split))
            theta[termID] <-
                attr(modelTerms, "start")[[i]](structure(theta[termID],
                                                         assign = split),
                                               gnmData)
        }
    }
    theta <- theta[uniq]

    if (attr(modelTerms, "intercept"))
        termAssign <- termAssign - 1

    prodList <- vector(mode = "list", length = nFactor)
    names(prodList) <- varLabels
    type <- type == "Special"

    varPredictors <- function(theta) {
        for (i in seqFactor) {
            prodList[[i]] <- .Call("submatprod", baseMatrix,
                                   theta[thetaID[[i]]],
                                   first[a[i]], nr, nc[i],
                                   PACKAGE = "gnm", NAOK = TRUE)
        }
        prodList
    }

    predictor <- function(varPredictors, term = FALSE) {
        if (term) {
            es <- lapply(attr(modelTerms, "predictor"), function(x) {
                do.call("bquote", list(x, gnmData))})
            tp <- matrix(sapply(es, eval, varPredictors), nr)
            colnames(tp) <- c("(Intercept)"[attr(modelTerms, "intercept")],
                              attr(modelTerms, "term.labels"))
            tp
        }
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
    tmpID <- tmpID[type[tmpID]]
    nCommon <- as.integer(nCommon[as.character(tmpID)])
    if (any(type))
        specialVarDerivs <- deriv(e, varLabels[type])
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
                if (type[fi]) {
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
                .Call("onecol", baseMatrix, as.double(unlist(v[fi])),
                      first[i1], lt[fi[1]], nr, as.integer(length(fi)),
                      PACKAGE = "gnm")
            }
        }
    }

    toolList <- list(start = theta, varPredictors = varPredictors,
                     predictor = predictor, localDesignFunction = localDesignFunction)
    if (x) toolList$termAssign <- termAssign[uniq]
    toolList
}
