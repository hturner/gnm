myBeliSNoInt <- function(x, peakage = ~1, pksharp = 1, left.S = 1,
                 constraint = c(min(x), max(x)), D = 5,
                 inst = NULL){
    call <- match.call()
    match <- match(c("peakage", "pksharp", "left.S"),
                   names(match.call()[-1]), 0)
    list(predictors = list(peakage = peakage,
         pksharp = substitute(pksharp),
         left.S = substitute(left.S)),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             SL <- paste(predLabels[3], sep = "")
             paste(" - exp(", predLabels[2], ") * ((",
                   SL, " * log(", SL, "/(",
                   SL, " + ", varLabels[1], " - ", predLabels[1], ")) + ",
                   varLabels[1], " - ", predLabels[1], ")/(",
                   SL, " * log(", SL, "/(", SL, " - ", D, ")) - ",
                   D, "))", sep = "")
         },
         call = as.expression(call),
         match = match
         )
}
class(myBeliSNoInt) <- "nonlin"
