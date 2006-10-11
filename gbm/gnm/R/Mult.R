Mult <- function(..., inst = NULL){
    if ("multiplicity" %in% names(match.call()[-1]))
        stop("multiplicity argument of Mult has been replaced by",
             "\"inst\" argument.")
    list(predictors = match.call(expand.dots = FALSE)[["..."]],
         term = function(predLabels, ...) {
             paste("(", paste(predLabels, collapse = ")*("), ")", sep = "")
         },
         call = as.expression(match.call()))
}
class(Mult) <- "nonlin"
