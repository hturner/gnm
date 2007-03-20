Bell <- function(x, slope = 1, endpoint = 1, side = "left",
                 constraint = ifelse(side == "right", max(x) + 1, min(x) - 1),
                 inst = NULL){
    call <- match.call()
    #call$slope <- call$endpoint <- NULL
    match <- match(c("slope", "endpoint"), names(match.call()[-1]), 0)
    list(predictors = list(slope = substitute(slope),
         endpoint = substitute(endpoint)),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             paste(predLabels[1], " * log(",
                   " -"[side == "right"], varLabels[1], " + ",
                   " -"[side == "left"], constraint,
                   " + exp(", predLabels[2], "))")
         },
         call = as.expression(call),
         match = match,
         start = function(theta){
             endpoint <- grep("[).]endpoint", names(theta))
             theta[endpoint] <- 0
             theta
         }
         )
}
class(Bell) <- "nonlin"
