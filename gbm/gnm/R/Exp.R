Exp <- function(multiplier){
    badCall <- charmatch(c("model.frame.default", "model.matrix.default"),
                         sapply(sys.calls(),
                                function(x) as.character(x[[1]])[1]))
    if (any(!is.na(badCall)))
        stop(paste("Exp terms are only valid in gnm models."))
    
    structure(deparse(substitute(multiplier)), class = "Exp")
}

