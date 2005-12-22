gnmTerms <- function(formula, eliminate)
{
    if (!is.null(substitute(e, list(e = eliminate)))) {
        tmp <- .Internal(update.formula(formula,
                                        substitute(~ -1 + I(e) + .,
                                                   list(e = eliminate))))
        environment(tmp) <- environment(formula)
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
                                  class = "Linear"))[any(
                                  c(intercept, !nonlinear))],
                   lapply(labelList[nonlinear],
                          function(term) eval(parse(text = term))))
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
