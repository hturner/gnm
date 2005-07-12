NonlinVariables <- function(functionCall) {
    varMethod <- paste(as.character(functionCall[[1]]), "Variables", sep = "")
    if (exists(varMethod)) {
        functionCall[[1]] <- as.name(varMethod)
        eval(functionCall)
    }
    else
        as.character(match.call(match.fun(functionCall[[1]]), functionCall,
                                expand.dots = FALSE)[["..."]])
}
