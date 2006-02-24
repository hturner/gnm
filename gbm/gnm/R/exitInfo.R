exitInfo <- function(object){
    conv <- object$converged
    if (conv)
        cat("Algorithm converged")
    else {
        cat("\nTolerance: ", object$tolerance, "\n")
        cat("\nAbsolute scores >= ",
            "tolerance * sqrt(tolerance + diag(information matrix)):\n\n")
        score <- abs(attr(conv, "score"))
        fail <- score >= attr(conv, "criterion")
        print(data.frame(abs.score = score,
                         criterion =  attr(conv, "criterion"))[fail,])
    }
}
