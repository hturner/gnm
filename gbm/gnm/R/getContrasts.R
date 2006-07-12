getContrasts <- function(model, sets = NULL,
                         nSets = if (is.list(sets)) length(sets) else 1,
                         dispersion = NULL,
                         use.eliminate = TRUE,
                         ...){
    result.as.list <- {
        (is.null(sets) && (nSets > 1)) ||
        (!is.null(sets) && is.list(sets))
    }
    coefs <- coef(model)
    l <- length(coefs)
    of.interest <- ofInterest(model)
    if (is.null(of.interest)) {
        if (l == model$eliminate) stop("model has no parameters of interest")
        of.interest <- (model$eliminate + 1):l
    }
    coefNames <- names(coefs)
    if (is.null(sets)) {
        of.interest <- ofInterest(model)
        if (is.null(of.interest)) {
            if (l == model$eliminate)
                stop("model has no parameters of interest")
            of.interest <- (model$eliminate + 1):l
        }
        sets <- relimp:::pickFrom(coefNames[of.interest], nSets, ...)
    }
    if (!is.list(sets)) sets <- list(sets)
    setLengths <- sapply(sets, length)
    if (all(setLengths == 0)) stop(
            "no non-empty parameter set specified")
    if (all(setLengths < 1.5)) stop(
            "for contrasts, at least 2 parameters are needed in a set")
    if (is.list(sets)) sets <- sets[setLengths > 1.5]
    if (any(setLengths < 1.5)) warning(
            "Sets with fewer than 2 parameters were dropped,")
    sets <- lapply(sets, function(x){
        if (is.numeric(x)) x <- coefNames[x]
        contr <- contr.sum(factor(x))
        rnames <- rownames(contr)
        contr <- rbind(contr[nrow(contr), ], contr[-nrow(contr), ])
        rownames(contr) <- rnames
        list(coefs = x, contr = contr)
    }
                   )
    coefMatrix <- lapply(sets, function(x){
        temp <- matrix(0, l, length(x$coefs))
        id <- match(x$coefs, coefNames)
        temp[id, 2:ncol(temp)] <- x$contr
        colnames(temp) <- x$coefs
        temp})
    Vcov <-  vcov(model, dispersion = dispersion,
                  use.eliminate = use.eliminate)
    result <- lapply(coefMatrix, function(x)
       {
        iden <- checkEstimable(model, x)
        if (any(!na.omit(iden))) {
            print(iden)
            cat("Note: not all of the specified contrasts in this set are",
                "estimable\n")
        }
        not.unestimable <- iden | is.na(iden)
        result <- se(model, x[, not.unestimable, drop = FALSE],
           checkEstimability = FALSE, Vcov = Vcov)
        relerrs <- NULL
        V <- NULL
        if (any(not.unestimable)){
            estimable.names <- names(not.unestimable)[not.unestimable]
            V <- Vcov[estimable.names, estimable.names, drop = FALSE]
        }
        if (sum(not.unestimable) > 2) {
            QVs <- qvcalc:::qvcalc(V)
            quasiSE <- sqrt(QVs$qvframe$quasiVar)
            result <- cbind(result, quasiSE)
            names(result)[1:2] <- c("estimate", "SE")
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
    if (result.as.list) result else result[[1]]
}
