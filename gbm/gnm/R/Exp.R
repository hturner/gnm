Exp <- function(factor){
    badCall <- charmatch(c("model.frame.default", "model.matrix.default"),
                         sapply(sys.calls(),
                                function(x) as.character(x[[1]])))
    if (!all(is.na(badCall))) {
        culprit <- gsub(".default", "",
                        as.character(sys.calls()
                                     [[min(badCall[!is.na(badCall)])]][[1]]))
        stop(paste(culprit,
                   "has called Exp() from the gnm package.\n", culprit,
                   "can only handle Exp terms",
                   "as part of the formula of a gnm object."))
    }
    
    structure(deparse(substitute(factor)), class = "Exp")
}

