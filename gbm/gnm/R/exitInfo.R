exitInfo <- function(object){
    conv <- object$converged
    if (conv)
        cat("Algorithm converged")
    else {
        cat("\nTolerance: ", object$tolerance, "\n")
        cat("\nAbsolute scores >= ",
            "tolerance * sqrt(tolerance + diag(information matrix)):\n")
        print(attr(conv, "score"))
    }
}
