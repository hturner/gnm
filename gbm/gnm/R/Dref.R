Dref <- function(..., delta = ~ 1){
    preds <- match.call(expand.dots = FALSE)[["..."]]
    n <- length(preds)
    preds <- c(delta = rep(list(delta), n), preds)
    common <- c(1:n, rep(n + 1, n))
    nf <- match(c("delta"), names(match.call()[-1]), 0)
    if ("formula" %in% names(match.call()[-1]))
        stop("formula argument of old plug-in has been renamed ",
             "\"delta\" in this function.")
    match <- c(rep(nf, n), 1:n)
    names(preds) <- c(rep("delta", n), rep("", n))

    list(predictors = preds,
         common = common,
         match = match,
         term = function(predLabels, ...){
             delta <- predLabels[1:n]
             gamma <- predLabels[-c(1:n)]
             paste("(((exp(", delta, "))/(",
                   paste("exp(", delta, ")", collapse = " + "),
                   "))*", gamma, ")", sep = "", collapse = " + ")},
         start = function(theta) {
             ifelse(attr(theta, "assign") == 3, 0.5, runif(length(theta)) - 0.5)
         },
         call = as.expression(match.call()))
}
class(Dref) <- "nonlin"
