Nonlin <- function(functionCall){
    badCall <- charmatch(c("model.frame.default", "model.matrix.default"),
                           sapply(sys.calls(),
                                  function(x) as.character(x[[1]])[1]))
    if (!all(is.na(badCall))) {
        culprit <- gsub(".default", "",
                        as.character(sys.calls()
                                     [[min(badCall[!is.na(badCall)])]][[1]]))
        stop(paste(culprit,
                   "has called Nonlin() from the gnm package.\n", culprit,
                   "can only handle Nonlin terms",
                   "as part of the formula of a gnm object."))
    }

    functionCall <- match.call()$functionCall
    envir <- environment(eval(functionCall[[1]]))
    varMethod <- paste(as.character(functionCall[[1]]), "Variables", sep = "")
    if (exists(varMethod, envir)) {
        varMethod <- as.call(c(as.name(varMethod), as.list(functionCall)[-1]))
        variables <- eval(substitute(varMethod), envir)
    }
    else
       variables <- as.character(match.call(match.fun(functionCall[[1]]),
                                            functionCall,
                                            expand.dots = FALSE)[["..."]])
    if (!length(variables))
        stop("No variables in term!")
    structure(variables, class = "Nonlin", call = functionCall)
}
