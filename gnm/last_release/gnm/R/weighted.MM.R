weighted.MM <- function(resp, conc){
    list(predictors = list(Vm = substitute(conc), K = 1),
         variables = list(substitute(resp), substitute(conc)),
         term = function(predLabels, varLabels) {
             pred <- paste("(", predLabels[1], "/(", predLabels[2],
                           " + ", varLabels[2], "))", sep = "")
             pred <- paste("(", varLabels[1], " - ", pred, ")/sqrt(",
                           pred, ")", sep = "")
         })
}
class(weighted.MM) <- "nonlin"
