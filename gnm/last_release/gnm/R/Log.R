Log <- function(expression, inst = NULL){
    list(predictors = list(substitute(expression)),
         term = function(predLabels, ...) {
             paste("log(", predLabels, ")", sep = "")
         },
         call = as.expression(match.call()),
         match = 1)
}
class(Log) <- "nonlin"
