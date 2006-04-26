gnmTerms <- function(formula, eliminate, data)
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
    fullTerms <- terms(formula, specials = "instances", keep.order = TRUE,
                       data = data)
    inst <- attr(fullTerms, "specials")$instances
    if (length(inst)) {
        instLabels <- rownames(attr(fullTerms, "factors"))[inst]
        instLabels2 <- sapply(instLabels, function(x) eval(parse(text = x)))
        formula <- deparse(formula, width.cutoff = 500)
        formula <- as.formula(mapply(gsub, instLabels, instLabels2,
                                     MoreArgs = list(formula, fixed = TRUE)))
        fullTerms <- terms(formula, keep.order = TRUE, data = data)
    }
    else
        attr(fullTerms, "specials") <- NULL
    
    if (is.empty.model(fullTerms))
        return(structure(formula, terms = fullTerms))

    labelList <- attr(fullTerms, "term.labels")
    intercept <- attr(fullTerms, "intercept")
    
    nonlinear <- is.element(seq(labelList),
                            grep("(Exp|Mult|Nonlin)[[:space:]]*\\(", labelList))
    labelList <- c(list(structure(c(intercept, labelList[!nonlinear]),
                                  class = "Linear", prefix = "", instance = "")
                        )[any(c(intercept, !nonlinear))],
                   lapply(labelList[nonlinear],
                          function(term) eval(parse(text = term))))
    prefixLabels <- sapply(labelList, attr, "prefix")
    instanceLabels <- sapply(labelList, attr, "instance")
    nonsense <- tapply(instanceLabels, prefixLabels, FUN = function(x)
                {nchar(x) && !identical(as.integer(x), seq(x))})
    if(any(nonsense))
        stop("Specified instances of ", prefixLabels[nonsense],
             " are not in sequence")
    constituentLabels <- sapply(labelList, attr, "constituentLabels")
    prefixLabels <- unlist(mapply(paste, paste(prefixLabels, instanceLabels,
                                               sep = ""),
                                  constituentLabels, sep = ""))

    labelList <- unlistOneLevel(labelList)
    offsetList <- lapply(labelList, attr, "offset")

    if (attr(fullTerms, "response") < 1) response <- NULL
    else response <- attr(fullTerms, "variables")[[2]]
    
    predictorOffset <- rownames(attr(fullTerms, "factors")
                                )[attr(fullTerms, "offset")]
    
    structure(reformulate(c(unlist(labelList), unlist(offsetList),
                            predictorOffset), response),
              terms = fullTerms,
              offset = offsetList,
              parsedLabels = labelList,
              prefixLabels = prefixLabels,
              .Environment = environment(formula))
}
