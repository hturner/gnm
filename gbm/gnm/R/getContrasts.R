getContrasts <- function(model, sets = NULL, nSets = 1, ...){
    if (is.null(sets)){
        if (!require(tcltk)) stop(
               "no parameter set specified, and tcltk not installed")
        if (!require(relimp)) stop(
               "the relimp package from CRAN needs to be installed")
        sets <- pickFrom(names(coef(model)), nSets,...)
        if (length(sets) == 0) stop("no parameter set(s) specified")
    }
    coefs <- coef(model)
    l <- length(coefs)
    if (!is.list(sets)) sets <- list(sets)
    nSets <- length(sets)
    sets <- lapply(sets, function(x){
        if (is.numeric(x)) x <- names(coefs)[x]
        list(coefs = x, contr = contr.sum(factor(x)))})
    coefMatrix <- lapply(sets, function(x){
        temp <- matrix(0, l, length(x$coefs))
        rownames(temp) <- names(coefs)
        temp[x$coefs, 1:(ncol(temp) - 1)] <- x$contr
        colnames(temp) <- x$coefs
        temp})
    lapply(coefMatrix, function(x){
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
        if (sum(not.unestimable) > 2 && require(qvcalc)) {
            not.unestimable <- names(not.unestimable)[not.unestimable]
            V <- vcov(model)[not.unestimable, not.unestimable, drop = FALSE]
            QVs <- qvcalc(V)
            quasi.se <- sqrt(QVs$qvframe$quasiVar)
            result <- cbind(result, quasi.se)
            relerrs <- round(100*c(min(QVs$relerrs), max(QVs$relerrs)), 1)
            relerrs <- paste(relerrs, "%", sep = "")
        }
        list(summary = result, relative.errors = relerrs)
    })
}
