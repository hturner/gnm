gnmConstrain <- function(model, constrain, start = coef(model), ...) {    
    call <- update(model, constrain = constrain, start = start,
                   evaluate = FALSE, ...)
    call$constrain <- do.call("substitute", list(constrain))
    call$start <- do.call("substitute", list(start))
    eval(call)
}
