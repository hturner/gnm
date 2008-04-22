Nonlin <- function(functionCall){
    checkCall()
    functionCall <- match.call()$functionCall
    if (deparse(functionCall[[1]]) %in% c("Dref", "MultHomog"))
        stop(deparse(functionCall[[1]]), " is now implemented as a ",
              "\"nonlin\" function")
    else
        warning("Plug-in functions are likely to be deprecated in future ",
                "versions. \nFunctions of class \"nonlin\" should be used ",
                "instead.", call. = FALSE)
         
    envir <- environment(eval(functionCall[[1]]))
    varMethod <- paste(functionCall[[1]], "Variables", sep = "")
    if (exists(varMethod, envir)) {
        varMethod <- as.call(c(as.name(varMethod), as.list(functionCall)[-1]))
        variables <- as.list(attr(terms(reformulate(eval(substitute(varMethod),
                                                         envir))),
                                  "variables"))[-1]
    }
    else
       variables <- match.call(match.fun(functionCall[[1]]), functionCall,
                               expand.dots = FALSE)[["..."]]
    if (!length(variables))
        stop("No variables in term!")

    Call <- deparse(sys.call())
        
    list(prefix = deparse(functionCall),
         matchID = 0,
         variables = variables,
         predvars = variables,
         varLabels =  Call,
         unitLabels = Call,
         block = 0,
         common = 1,
         classID = "Nonlin",
         NonlinID = "Nonlin",
         predictor = paste("`", Call, "`", sep = ""))
}
class(Nonlin) <- "nonlin"
