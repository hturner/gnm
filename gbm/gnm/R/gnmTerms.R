gnmTerms <- function(formula, eliminate)
{
    if (!is.null(eliminate)) {
        tmp <- .Internal(update.formula(formula,
                                        substitute(~ -1 + e + .,
                                                   list(e = eliminate[[2]]))))
        formula <- formula(terms.formula(tmp, simplify = TRUE,
                                         keep.order = TRUE))
    }
    fullTerms <- terms(formula, keep.order = TRUE)
    if (is.empty.model(fullTerms))
        return(structure(formula, terms = fullTerms))

    labelList <- attr(fullTerms, "term.labels")
    intercept <- attr(fullTerms, "intercept")
    nonlinear <- is.element(seq(labelList),
                            grep("(Mult|Nonlin)[[:space:]]*\\(", labelList))
    labelList <- c(list(structure(c(intercept, labelList[!nonlinear]),
                                  class = "Linear"))[any(!nonlinear|intercept)],
                   lapply(labelList[nonlinear],
                          function(term) eval(parse(text = term))))
    labelList <- prefixList <- unlistOneLevel(labelList)
    
    classIndex <- sapply(labelList, class)
    multNo <- cumsum(classIndex == "Mult")
    for (i in seq(labelList))
        prefixList[[i]] <-
            switch(classIndex[[i]],
                   "Mult" = paste("Mult", multNo[i], ".Factor",
                   seq(labelList[[i]]), ".", sep = ""),
                   "Nonlin" = deparse(attr(labelList[[i]], "call")),
                   "")

    labelList <- unlistOneLevel(labelList)
    offsetList <- lapply(labelList, attr, "offset")

    extraData <- lapply(labelList, attr, "extraData")
    extraData <- do.call("cbind", extraData[!sapply(extraData, is.null)])

    if (attr(fullTerms, "response") < 1) response <- ""
    else response <- evalq(attr(fullTerms, "variables")[[2]])
    
    predictorOffset <- sapply((attr(fullTerms,
                                    "variables")[-1])[attr(fullTerms,
                                                           "offset")], deparse)
    
    structure(reformulate(c(unlist(labelList), unlist(offsetList),
                            predictorOffset, unlist(colnames(extraData))),
                          response),
              terms = fullTerms,
              offset = offsetList,
              parsed.labels = labelList,
              prefix.labels = unlist(prefixList),
              response = attr(fullTerms, "response"),
              extraData = extraData)
}
