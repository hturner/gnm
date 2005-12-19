getContrasts <- function(model, sets = NULL, nSets = 1, ...){
    if (is.null(sets)){
        if (!require(tcltk)) stop(
               "no parameter set specified, and tcltk not installed")
        if (!require(relimp)) stop(
               "the relimp package from CRAN needs to be installed")
        sets <- pickFrom(names(coef(model)), nSets,...)
    }
    setLengths <- sapply(sets, length)
    if (all(setLengths == 0)) stop(
            "no non-empty parameter set specified")
    if (all(setLengths < 1.5)) stop(
            "for contrasts, at least 2 parameters are needed in a set")
    sets <- sets[setLengths > 1.5]
    if (any(setLengths < 1.5)) warning(
            "Sets with fewer than 2 parameters were dropped,")
    coefs <- coef(model)
    l <- length(coefs)
    if (!is.list(sets)) sets <- list(sets)
    nSets <- length(sets)
    sets <- lapply(sets, function(x){
        if (is.numeric(x)) x <- names(coefs)[x]
        contr <- contr.sum(factor(x))
        rnames <- rownames(contr)
        contr <- rbind(contr[nrow(contr), ], contr[-nrow(contr), ])
        rownames(contr) <- rnames
        list(coefs = x, contr = contr)
    }
                   )
    coefMatrix <- lapply(sets, function(x){
        temp <- matrix(0, l, length(x$coefs))
        rownames(temp) <- names(coefs)
        temp[x$coefs, 2:ncol(temp)] <- x$contr
        colnames(temp) <- x$coefs
        temp})
    lapply(coefMatrix, function(x)
       {
        iden <- checkEstimable(model, x)
        if (any(!na.omit(iden))) {
            print(iden)
            cat("Note: not all of the specified contrasts in this set are",
                "estimable\n")
        }
        not.unestimable <- iden | is.na(iden)
        result <- se(model, x[, not.unestimable, drop = FALSE],
           checkEstimability = FALSE)
        relerrs <- NULL
        V <- NULL
        if (any(not.unestimable)){
            estimable.names <- names(not.unestimable)[not.unestimable]
            V <- vcov(model)[estimable.names, estimable.names, drop = FALSE]
        }
        if (sum(not.unestimable) > 2 && require(qvcalc)) {
            QVs <- qvcalc(V)
            quasiSE <- sqrt(QVs$qvframe$quasiVar)
            result <- cbind(result, quasiSE)
            result$quasiVar <- QVs$qvframe$quasiVar
            relerrs <- QVs$relerrs

        }
        return(structure(list(covmat = V,
                              qvframe = result,
                              relerrs = relerrs,
                              modelcall = model$call),
                         class = "qv")
               )
    }
           )
}
