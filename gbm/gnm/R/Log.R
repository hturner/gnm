Log <- function(expression, inst = NULL){
    list(predictors = list(substitute(expression)),
         term = function(predLabels, ...) {
             paste("log(", predLabels, ")", sep = "")
         },
         call = as.expression(match.call()))
}
class(Log) <- "nonlin"
