Plogis <- Expit <- function(x, inst = NULL){
    list(variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             paste("1/(1 + exp(-", varLabels, "))", sep = "")
         },
         call = as.expression(match.call()),
         match = 1)
}
class(Plogis) <- class(Expit) <- "nonlin"
