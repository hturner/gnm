Mult <- function(..., inst = 1){
    badCall <- charmatch(c("model.frame.default", "model.matrix.default"),
                           sapply(sys.calls(),
                                  function(x) as.character(x[[1]])[1]))
    if (any(!is.na(badCall)))
        stop(paste("Mult terms are only valid in gnm models."))
    
    factorList <- as.list(as.character((match.call(expand.dots = FALSE))[[2]]))
    for (i in grep("Exp[[:space:]]*\\(", unlist(factorList)))
        factorList[[i]] <- eval(parse(text = factorList[[i]]))
    factorList <- lapply(factorList, function(x) {
        xTerms <- terms(as.formula(paste("~", x)))
        if (!is.null(attr(xTerms, "offset"))) {
            structure(paste(c(attr(xTerms, "term.labels"),
                              attr(xTerms, "intercept")), collapse = " + "),
                      offset = paste(sapply(lapply((attr(xTerms,
                      "variables")[-1])[attr(xTerms, "offset")], "[[", 2),
                      deparse), collapse = " + "))
        }
        else
            x
    })

    Call <- sys.call()
    Call$inst <- NULL
    
    structure(factorList, class = "Mult",
              prefix = deparse(Call)[1], instance = paste("", inst, sep = ""),
              constituentLabels = paste(".", seq(factorList), sep = ""))
}

