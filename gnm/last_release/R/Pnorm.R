Pnorm <- Phi <- function(x, inst = NULL){
    list(variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             paste("pnorm(", varLabels, ")", sep = "")
         },
         call = as.expression(match.call()),
         match = 1)
}
class(Pnorm) <- class(Phi) <- "nonlin"
