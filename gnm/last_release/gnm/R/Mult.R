Mult <- function(..., inst = NULL){
    if ("multiplicity" %in% names(match.call()[-1]))
        stop("multiplicity argument of Mult has been replaced by",
             "\"inst\" argument.")
    dots <- match.call(expand.dots = FALSE)[["..."]]
    list(predictors = dots,
         term = function(predLabels, ...) {
             paste("(", paste(predLabels, collapse = ")*("), ")", sep = "")
         },
         call = as.expression(match.call()),
         match = seq(dots))
}
class(Mult) <- "nonlin"
