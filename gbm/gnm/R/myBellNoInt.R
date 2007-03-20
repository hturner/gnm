myBellNoInt <- function(x, peakage = 1, pksharp = 1, left.ep = 1,
                 right.ep = 1, constraint = c(min(x) - 1, max(x) + 1),
                 inst = NULL){
    call <- match.call()
    match <- match(c("peakage", "pksharp", "left.ep", "right.ep"),
                   names(match.call()[-1]), 0)
    list(predictors = list(peakage = substitute(peakage),
         pksharp = substitute(pksharp),
         left.ep = substitute(left.ep), right.ep = substitute(right.ep)),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             SL <- paste("(", predLabels[1], " - ", constraint[1], " + 1e-5 + exp(",
                         predLabels[3], "))", sep = "")
             SR <- paste("(", constraint[2], " + 1e-5 + exp(", predLabels[4], ") - ",
                         predLabels[1], ")", sep = "")
             paste(" - exp(", predLabels[2], ") * ((",
                   SL, " * log(", SL, "/(",
                   SL, " + ", varLabels[1], " - ", predLabels[1], ")) + ",
                   SR, " * log(", SR, "/(",
                   SR, " - ", varLabels[1], " + ", predLabels[1], ")))/(",
                   SL, " * log(", SL, "/(", SL, " - 5)) + ",
                   SR, " * log(", SR, "/(", SR, " + 5))))", sep = "")
         },
         call = as.expression(call),
         match = match
         )
}
class(myBellNoInt) <- "nonlin"
