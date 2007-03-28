myBeli <- function(x, peakage = ~1, peak.ht = 1, pksharp = 1, left.ep = 1,
                 constraint = c(min(x), max(x)), D = 5,
                 inst = NULL){
    call <- match.call()
    match <- match(c("peakage", "peak.ht", "pksharp", "left.ep"),
                   names(match.call()[-1]), 0)
    list(predictors = list(peakage = peakage,
         peak.ht = substitute(peak.ht), pksharp = substitute(pksharp),
         left.ep = substitute(left.ep)),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             SL <- paste("(", predLabels[1], " - ", constraint[1],
                         " + 1e-5 + 1e5/(1 + exp(-",
                         predLabels[4], ")))", sep = "")
             paste(predLabels[2], " - exp(", predLabels[3], ") * ((",
                   SL, " * log(", SL, "/(",
                   SL, " + ", varLabels[1], " - (", predLabels[1], "))) + ",
                   varLabels[1], " - (", predLabels[1], "))/(",
                   SL, " * log(", SL, "/(", SL, " - ", D, ")) - ",
                   D, "))", sep = "")
         },
         call = as.expression(call),
         match = match
         )
}
class(myBeli) <- "nonlin"
