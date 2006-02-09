model.frame.gnm <- function (formula, ...) 
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 
        0)]
    if (length(nargs) || is.null(formula$model)) {
        args <- formula$call
        args$method <- "model.frame"
        args[names(nargs)] <- nargs
        args[[1]] <- as.name("list")
        env <- environment(formula$terms)
        if (is.null(env)) 
            env <- parent.frame()
        do.call("gnm", eval(args, env))
    }
    else formula$model
}
        
        
    
