Mult <- function(..., multiplicity = 1){
    badCall <- charmatch(c("model.frame.default", "model.matrix.default"),
                           sapply(sys.calls(),
                                  function(x) as.character(x[[1]])[1]))
    if (!all(is.na(badCall))) {
        culprit <- gsub(".default", "",
                        as.character(sys.calls()
                                     [[min(badCall[!is.na(badCall)])]][[1]]))
        stop(paste(culprit,
                   "has called Mult() from the gnm package.\n", culprit,
                   "can only handle Mult terms",
                   "as part of the formula of a gnm object."))
    }
    
    factorList <- as.list(as.character((match.call(expand.dots = FALSE))[[2]]))
    for (i in grep("Exp[[:space:]]*\\(", unlist(factorList)))
        factorList[[i]] <- eval(parse(text = factorList[[i]]))
    factorList <- lapply(factorList, function(x) {
        xTerms <- terms(as.formula(paste("~", x)))
        if (!is.null(attr(xTerms, "offset"))) {
            structure(paste(attr(xTerms, "term.labels"), collapse = " + "),
                      offset = paste(sapply(lapply((attr(xTerms,
                      "variables")[-1])[attr(xTerms, "offset")], "[[", 2),
                      deparse), collapse = " + "))
        }
        else
            x
    })
    class(factorList) <- "Mult"
    rep(list(factorList), multiplicity)
}

