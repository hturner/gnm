Nonlin <- function(functionCall){
    .Defunct(msg = paste("'Nonlin' is defunct.",
             "\nUse functions of class \"nonlin\" instead.",
             "\nSee ?nonlin.function for more details."))
}
class(Nonlin) <- "nonlin"
