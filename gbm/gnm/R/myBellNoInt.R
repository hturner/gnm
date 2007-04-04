myBellNoInt <- function(x, peakage = ~1, pksharp = 1, left.ep = 1,
                 right.ep = 1, constraint = c(min(x), max(x)),
                 inst = NULL, D = 5){
    call <- match.call()
    match <- match(c("peakage", "pksharp", "left.ep", "right.ep"),
                   names(match.call()[-1]), 0)
    peakage.start <- mean(x) #location of peak
    ## pksharp.start <- log(0.5) #to give pksharp = 0.5, when = D = 5
    pksharp.start <- 0
    ep.start <- -log(1e5/(1 - 1e-5) - 1) #to give plus 1 year
    list(predictors = list(peakage = peakage,
         pksharp = substitute(pksharp),
         left.ep = substitute(left.ep), right.ep = substitute(right.ep)),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             SL <- paste("(", predLabels[1], " - ", constraint[1],
                         " + 1e-5 + 1e5/(1 + exp(-",
                         predLabels[3], ")))", sep = "")
             SR <- paste("(", constraint[2],
                         " + 1e-5 + 1e5/(1 + exp(-", predLabels[4], ")) - ",
                         predLabels[1], ")", sep = "")
             paste(" - exp(", predLabels[2], ") * ((",
                   SL, " * log(", SL, "/(",
                   SL, " + ", varLabels[1], " - (", predLabels[1], "))) + ",
                   SR, " * log(", SR, "/(",
                   SR, " - ", varLabels[1], " + (", predLabels[1], "))))/(",
                   SL, " * log(", SL, "/(", SL, " - ", D, ")) + ",
                   SR, " * log(", SR, "/(", SR, " + ", D, "))))", sep = "")
         },
         call = as.expression(call),
         match = match,
         start = function(theta) {
             nm <- names(theta)
             gnmData <- get("gnmData", parent.frame())
             formula <- update(peakage, y ~ .)
             peakage.start <- lm(formula,
                                 data = cbind(gnmData, y = peakage.start))
             peakage <- grep("[)]\\.?peakage", nm)
             pksharp <- grep("[)]\\.?pksharp", nm)
             ep <- grep("[.]ep", nm)
             theta[peakage] <- coef(peakage.start)
             theta[pksharp] <- pksharp.start
             theta[ep] <- ep.start
             theta
         }
         )
}
class(myBellNoInt) <- "nonlin"
