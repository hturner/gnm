Nonlin <- function(functionCall){
    .Defunct(msg = paste("'Nonlin' is defunct.",
             "\nUse functions of class \"nonlin\" instead.",
             "\nSee ?nonlin.function for more details."))
}
class(Nonlin) <- "nonlin"

getModelFrame <- function() {
    .Defunct(msg = paste("'getModelFrame' is deprecated as it was designed to ",
             "work with the old plug-in architecture for nonlinear terms."))
}

qrSolve <- function(A, b, rank = NULL, ...) {
    .Defunct(msg = paste("'qrSolve' is deprecated as it is no longer used ",
             "by gnm."))
}
