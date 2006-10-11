Inv <- function(expression, inst = NULL){
    list(predictors = list(substitute(expression)),
         term = function(predLabels, ...) {
             paste("(", predLabels, ")^-1", sep = "")
         },
         call = as.expression(match.call()))
}
class(Inv) <- "nonlin"
