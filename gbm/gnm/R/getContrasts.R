getContrasts <- function(model, sets = NULL, nsets = 1, ...){
    if (is.null(sets)){
        sets <- pickFrom(names(coef(model)), nsets,...)
        if (length(sets) == 0) stop("no parameter set(s) specified")
    }
    coefs <- coef(model)
    l <- length(coefs)
    if (!is.list(sets)) sets <- list(sets)
    nsets <- length(sets)
    sets <- lapply(sets, function(x){
        if (is.numeric(x)) x <- names(coefs)[x]
        list(coefs = x, contr = contr.sum(factor(x)))})
    cmatrix <- lapply(sets, function(x){
        temp <- matrix(0, l, length(x$coefs))
        rownames(temp) <- names(coefs)
        temp[x$coefs, 1:(ncol(temp) - 1)] <- x$contr
        colnames(temp) <- x$coefs
        temp})
    lapply(cmatrix, function(x){
        iden <- checkEstimable(model, x)
        if (any(!na.omit(iden))) {
            print(iden)
            cat("Note: not all of the specified contrasts in this set are",
                "estimable\n")
        }
        not.unestimable <- iden | is.na(iden)
        se(model, x[, not.unestimable, drop = FALSE],
           check.estimability = FALSE)
    })
}
