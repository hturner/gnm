Logistic <- function(x, inst = NULL){
    list(predictors = list(Asym = 1, xmid = 1, scal = 1),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             paste(predLabels[1], "/(1 + exp((", predLabels[2], "-",
                   varLabels[1], ")/", predLabels[3], "))")
         },
         start = function(theta){
             c(NA, mean(x), sd(x))
         }
         )
}
class(Logistic) <- "nonlin"
