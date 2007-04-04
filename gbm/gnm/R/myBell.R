myBell <- function(x, peakage = 1, peak.ht = 1, pksharp = 1, left.ep = 1,
                 right.ep = 1, constraint = c(min(x), max(x)), D = 5,
                 inst = NULL){
    call <- match.call()
    match <- match(c("peakage", "peak.ht", "pksharp", "left.ep", "right.ep"),
                   names(match.call()[-1]), 0)
    list(predictors = list(peakage = substitute(peakage),
         peak.ht = substitute(peak.ht), pksharp = substitute(pksharp),
         left.ep = substitute(left.ep), right.ep = substitute(right.ep)),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             SL <- paste("(", predLabels[1], " - ", constraint[1],
                         " + 1e-5 + 1e5/(1 + exp(-",
                         predLabels[4], ")))", sep = "")
             SR <- paste("(", constraint[2],
                         " + 1e-5 + 1e5/(1 + exp(-", predLabels[5], ")) - ",
                         predLabels[1], ")", sep = "")
             paste(predLabels[2], " - exp(", predLabels[3], ") * ((",
                   SL, " * log(", SL, "/(",
                   SL, " + ", varLabels[1], " - (", predLabels[1], "))) + ",
                   SR, " * log(", SR, "/(",
                   SR, " - ", varLabels[1], " + (", predLabels[1], "))))/(",
                   SL, " * log(", SL, "/(", SL, " - ", D, ")) + ",
                   SR, " * log(", SR, "/(", SR, " + ", D, "))))", sep = "")
         },
         call = as.expression(call),
         match = match
         )
}
class(myBell) <- "nonlin"
