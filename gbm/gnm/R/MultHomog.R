MultHomog <- function(..., inst = NULL){
    dots <- match.call(expand.dots = FALSE)[["..."]]
    list(predictors = dots,
         common = rep(1, length(dots)),
         term = function(predLabels, ...) {
             paste("(", paste(predLabels, collapse = ")*("), ")", sep = "")
         },
         call = as.expression(match.call()),
         match = rep(0, length(dots)))
}
class(MultHomog) <- "nonlin"
