Logit <- function(expression, inst = NULL){
    list(predictors = list(substitute(expression)),
         term = function(predLabels, ...) {
             paste("log((", predLabels, ")/(1 - (", predLabels, ")))", sep = "")
         },
         call = as.expression(match.call()))
}
class(Logit) <- "nonlin"
