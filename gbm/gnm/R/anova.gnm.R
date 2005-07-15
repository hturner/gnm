anova.gnm <- function (object, ..., dispersion = NULL, test = NULL) 
{
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) 
        rep(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named)) 
        warning("the following arguments to 'anova.gnm' are invalid and dropped: ", 
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.gnm <- unlist(lapply(dotargs, function(x) inherits(x, "gnm")))
    dotargs <- dotargs[is.gnm]
    if (length(dotargs) > 0) 
        return(anova.glmlist(c(list(object), dotargs), dispersion = dispersion, 
            test = test))
    varlist <- attr(object$terms, "variables")
    x <- model.matrix(object)
    varseq <- attr(x, "assign")
    nvars <- max(0, varseq)
    resdev <- resdf <- NULL
    if (nvars > 0) {
        for (i in 0:(nvars - 1)) {
            constrain <- replace(object$constrain, varseq > i, TRUE)
            fit <- update(object, constrain = constrain, verbose = FALSE)
            resdev <- c(resdev, fit$deviance)
            resdf <- c(resdf, fit$df.residual)
        }
        resdf <- c(resdf, object$df.residual)
        resdev <- c(resdev, object$deviance)
        table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))), 
        resdf, resdev)
    }
    else
        table <- data.frame(NA, NA, object$df.residual, object$deviance)
    dimnames(table) <- list(c("NULL", attr(object$terms, "term.labels")),
                            c("Df", "Deviance", "Resid. Df", "Resid. Dev"))
    title <- paste("Analysis of Deviance Table", "\n\nModel: ", 
                   object$family$family, ", link: ", object$family$link, 
                   "\n\nResponse: ", as.character(varlist[-1])[1],
                   "\n\nTerms added sequentially (first to last)\n\n", 
                   sep = "")
    df.dispersion <- Inf
    if (is.null(dispersion)) {
        dispersion <- summary(object, dispersion = dispersion)$dispersion
        df.dispersion <- if (dispersion == 1) 
            Inf
        else object$df.residual
    }
    if (!is.null(test)) 
        table <- stat.anova(table = table, test = test, scale = dispersion, 
            df.scale = df.dispersion, n = NROW(x))
    structure(table, heading = title, class = c("anova", "data.frame"))
}
