getContrasts <- function(model, sets = NULL, nsets = 1, ...){
    if (is.null(sets)){
        require(tcltk)
        require(relimp)
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
        temp[x$coefs, 1:(ncol(temp)-1)] <- x$contr
        colnames(temp) <- x$coefs
        temp})
    all.contrasts <- matrix(NA, l, 0)
    for (cmat in cmatrix){
        all.contrasts <- cbind(all.contrasts, cmat)
    }
    iden <- checkIdentifiability(model, all.contrasts)
    if (any(!iden)){
        print(iden)
        cat("Error: some contrasts not identified in the model\n")
        return(NULL)
    }
    lapply(cmatrix, function(x){
        se(model, x, check.identifiability = FALSE)
    })
}
