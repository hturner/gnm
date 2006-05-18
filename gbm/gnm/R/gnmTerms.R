gnmTerms <- function(formula, eliminate = NULL, data = NULL)
{
    if (!is.null(substitute(e, list(e = eliminate)))) {
        tmp <- .Internal(update.formula(formula,
                                        substitute(~ -1 + e + .,
                                                   list(e = eliminate))))
        environment(tmp) <- environment(formula)
        data <- data[!names(data) %in% deparse(eliminate)]
        formula <- formula(terms.formula(tmp, simplify = TRUE,
                                         keep.order = TRUE, data = data))
    }
    fullTerms <- terms(formula, keep.order = TRUE, data = data)
    if (is.empty.model(fullTerms))
        return(structure(formula, terms = fullTerms))

    labelList <- as.list(attr(fullTerms, "term.labels"))
    if (attr(fullTerms, "intercept") == 1)
        labelList <- c("1", labelList)

    for (i in seq(labelList)) {
        if (length(grep("(Mult|Nonlin)[[:space:]]*\\(", labelList[[i]])))
            labelList[[i]] <- eval(parse(text = labelList[[i]]))
        else
            labelList[[i]] <- list(structure(labelList[[i]], class = "Linear"))
    }

    multiplicity <- sapply(labelList, function(x)
                           ifelse(is.list(x), length(x), 1))
    if (any(multiplicity > 1))
        termsID <- rep(seq(multiplicity), multiplicity)
    else
        termsID <- NULL
    labelList <- prefixList <- unlistOneLevel(labelList)
    
    classIndex <- sapply(labelList, class)
    multNo <- cumsum(classIndex == "Mult")
    for (i in seq(labelList))
        prefixList[[i]] <-
            switch(classIndex[[i]],
                   "Mult" = paste("Mult", multNo[i], ".Factor",
                   seq(labelList[[i]]), ".", sep = ""),
                   "Nonlin" = deparse(attr(labelList[[i]], "call"))[1],
                   "")

    labelList <- unlistOneLevel(labelList)
    offsetList <- lapply(labelList, attr, "offset")

    if (attr(fullTerms, "response") < 1) response <- NULL
    else response <- evalq(attr(fullTerms, "variables")[[2]])
    
    predictorOffset <- sapply((attr(fullTerms,
                                    "variables")[-1])[attr(fullTerms,
                                                           "offset")], deparse)
    
    structure(reformulate(c(unlist(labelList), unlist(offsetList),
                            predictorOffset), response),
              terms = fullTerms,
              termsID = termsID,
              offset = offsetList,
              parsedLabels = labelList,
              prefixLabels = unlist(prefixList),
              .Environment = environment(formula))
}
