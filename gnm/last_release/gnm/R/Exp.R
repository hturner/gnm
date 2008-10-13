Exp <- function(expression, inst = NULL){
    list(predictors = list(substitute(expression)),
         term = function(predLabels, ...) {
             paste("exp(", predLabels, ")", sep = "")
         },
         call = as.expression(match.call()),
         match = 1)
}
class(Exp) <- "nonlin"
