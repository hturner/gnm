myBeliNoInt <- function(x, peakage = ~1, pksharp = 1, left.ep = 1,
                 constraint = c(min(x), max(x)), D = 5,
                 inst = NULL){
    call <- match.call()
    match <- match(c("peakage", "pksharp", "left.ep"),
                   names(match.call()[-1]), 0)
    peakage.start <- mean(x) #location of peak
    pksharp.start <- log(0.5) #to give pksharp = 0.5
    left.ep.start <- -log(1e5/(1 - 1e-5) - 1) #to give plus 1 year
    list(predictors = list(peakage = peakage,
         pksharp = substitute(pksharp),
         left.ep = substitute(left.ep)),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             SL <- paste("(", predLabels[1], " - ", constraint[1],
                         " + 1e-5 + 1e5/(1 + exp(-",
                         predLabels[3], ")))", sep = "")
             paste(" - exp(", predLabels[2], ") * ((",
                   SL, " * log(", SL, "/(",
                   SL, " + ", varLabels[1], " - (", predLabels[1], "))) + ",
                   varLabels[1], " - (", predLabels[1], "))/(",
                   SL, " * log(", SL, "/(", SL, " - ", D, ")) - ",
                   D, "))", sep = "")
         },
         call = as.expression(call),
         match = match,
         start = function(theta) {
             nm <- names(theta)
             gnmData <- get("gnmData", parent.frame())
             formula <- update(peakage, y ~ .)
             peakage.start <- lm(formula,
                                 data = cbind(gnmData, y = peakage.start))
             peakage <- grep("\)\.?peakage", nm)
             pksharp <- grep("\)\.?pksharp", nm)
             left.ep <- grep("\)\.?left.ep", nm)
             theta[peakage] <- coef(peakage.start)
             theta[pksharp] <- pksharp.start
             theta[left.ep] <- left.ep.start
             theta
         }
         )
}
class(myBeliNoInt) <- "nonlin"
