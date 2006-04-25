Exp <- function(expression, inst = NULL){
    badCall <- charmatch(c("model.frame.default", "model.matrix.default"),
                         sapply(sys.calls(),
                                function(x) as.character(x[[1]])[1]))
    if (!all(is.na(badCall)))
        stop(paste("Exp terms are only valid in gnm models."))

    Call <- sys.call()
    Call$inst <- NULL
    
    structure(deparse(substitute(expression)), class = "Exp",
              prefix = deparse(Call)[1], instance = paste("", inst, sep = ""))
}

