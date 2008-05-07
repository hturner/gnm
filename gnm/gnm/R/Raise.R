Raise <- function(expression, power = 1, inst = NULL){
    list(predictors = list(substitute(expression)),
         term = function(predLabels, ...) {
             paste("(", predLabels, ")^", power, sep = "")
         },
         call = as.expression(match.call()),
         match = 1)
}
class(Raise) <- "nonlin"
