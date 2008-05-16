getContrasts <- function(model, set = NULL,
                         refLevel = "first",
                         dispersion = NULL,
                         use.eliminate = TRUE,
                         ...){
    coefs <- coef(model)
    l <- length(coefs)
    of.interest <- ofInterest(model)
    if (is.null(of.interest)) {
        if (l == model$eliminate) stop("Model has no parameters of interest")
        of.interest <- (model$eliminate + 1):l
    }
    coefNames <- names(coefs)
    if (is.null(set)) {
        of.interest <- ofInterest(model)
        if (is.null(of.interest)) {
            if (l == model$eliminate)
                stop("Model has no parameters of interest")
            of.interest <- (model$eliminate + 1):l
        }
        set <- unlist(relimp:::pickFrom(coefNames[of.interest], 1, ...))
    }
    setLength <- length(set)
    if (setLength == 0) stop(
            "No non-empty parameter set specified")
    if (setLength < 1.5) stop(
            "For contrasts, at least 2 parameters are needed in a set")
    if (is.numeric(set)) set <- coefNames[set]

    if (is.numeric(refLevel)){
        if (length(refLevel) == 0) stop("The specified refLevel has zero length")
        if (length(refLevel) == 1){
            if (refLevel %in% seq(setLength)) {
                temp <- rep(0, setLength)
                temp[refLevel] <- 1
                refLevel <- temp
            } else stop("The specified refLevel is out of range")
        }
        if (length(refLevel) != setLength) stop(
                  "The specified refLevel has the wrong length")
        if ((sum(refLevel) - 1) ^ 2 > 1e-10) stop(
                                   "The refLevel weights do not sum to 1")
    }
    if (identical(refLevel, "first")) refLevel <- c(1, rep(0, setLength - 1)) else
    if (identical(refLevel, "last")) refLevel <- c(rep(0, setLength - 1), 1) else
    if (identical(refLevel,"mean")) refLevel <- rep(1/setLength, setLength)

    contr <- diag(rep(1, setLength))
    contr <- contr - refLevel
    rownames(contr) <- set
    set <- list(coefNames = set, contr = contr)

    coefMatrix <- matrix(0, l, length(set$coefNames))
    id <- match(set$coefNames, coefNames)

    coefMatrix[id, ] <- set$contr
    colnames(coefMatrix) <- set$coefNames

    Vcov <-  vcov(model, dispersion = dispersion,
                  use.eliminate = use.eliminate)

    iden <- checkEstimable(model, coefMatrix)
    if (any(!na.omit(iden))) {
        cat("Estimability of the specified contrasts:\n")
        print(iden)
        if (all(!na.omit(iden))) stop(
                 "None of the specified contrasts is estimable")
        cat("Note: not all of the specified contrasts in this set are",
            "estimable\n")
    }
    not.unestimable <- iden | is.na(iden)
    result <- se(model, coefMatrix[, not.unestimable, drop = FALSE],
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
