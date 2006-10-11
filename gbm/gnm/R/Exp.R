Exp <- function(expression, inst = NULL){
    list(predictors = list(substitute(expression)),
         term = function(predLabels, ...) {
             paste("exp(", predLabels, ")", sep = "")
         },
         call = as.expression(match.call()))
}
class(Exp) <- "nonlin"
