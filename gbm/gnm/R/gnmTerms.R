gnmTerms <- function(formula)
{
    fullTerms <- terms(formula, keep.order = TRUE)
    if (is.empty.model(fullTerms))
        return(structure(formula, terms = fullTerms))

    labelList <- as.list(attr(fullTerms, "term.labels"))
    if (attr(fullTerms, "intercept") == 1) labelList <- c("1", labelList)
    for (i in grep("(Dref|Mult|Nonlin)[[:space:]]*\\(", unlist(labelList)))
        labelList[[i]] <- eval(parse(text = labelList[[i]]))
    labelList <- prefixList <- unlistOneLevel(labelList)
    
    classIndex <- sapply(labelList, class)
    multNo <- cumsum(classIndex == "Mult")
    for (i in seq(labelList))
        prefixList[[i]] <-
            switch(classIndex[[i]],
                   "Dref" = paste("Dref(", paste(labelList[[i]],
                   collapse = ","),")", sep = ""),
                   "Mult" = paste("Mult", multNo[i], ".Factor",
                   seq(labelList[[i]]), ".", sep = ""),
                   "Nonlin" = deparse(attr(labelList[[i]], "call")),
                   "")

    labelList <- unlistOneLevel(labelList)
    offsetList <- lapply(labelList, attr, "offset")

    extraData <- lapply(labelList, attr, "extraData")
    extraData <- do.call("cbind", extraData[!sapply(extraData, is.null)])

    response <- ifelse(attr(fullTerms, "response") < 1, "",
                             as.character(attr(fullTerms, "variables")[2]))
    predictorOffset <- sapply((attr(fullTerms,
                                    "variables")[-1])[attr(fullTerms,
                                                           "offset")], deparse)
    
    structure(as.formula(paste(response, "~",
                               paste(c(unlist(labelList), unlist(offsetList),
                                       predictorOffset,
                                       unlist(colnames(extraData))),
                                       collapse = "+"))),
              terms = fullTerms,
              offset = offsetList,
              parsed.labels = labelList,
              prefix.labels = unlist(prefixList),
              response = attr(fullTerms, "response"),
              extraData = extraData)
}
