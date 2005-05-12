Nonlin <- function(functionCall, data = NULL){
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
    variables <- match.call(match.fun(functionCall[[1]]), functionCall,
                            expand.dots = FALSE)[["..."]]
    structure(as.character(variables), class = "Nonlin", call = functionCall,
              extraData = substitute(data))
}
