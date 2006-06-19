se <- function(model, estimate = ofInterest(model),
               checkEstimability = TRUE, Vcov = NULL, dispersion = NULL, ...){
    if (!inherits(model, "gnm")) stop("model is not of class \"gnm\"")
    if (!is.null(Vcov) && !is.null(dispersion)) {
        Vcov <- Vcov * dispersion
    }
    coefs <- coef(model)
    l <- length(coefs)
    eliminate <- model$eliminate
    coefNames <- names(coefs)
    if (identical(estimate, "[?]"))
        estimate <- pickCoef(model, subset = (eliminate + 1):l,
                             title = paste("Estimate standard errors",
                             "for one or more gnm coefficients"))
    if (is.null(estimate))
        return(data.frame(coef(summary(model)))[, 1:2])
    else {
        if (is.character(estimate))
            estimate <- match(estimate, coefNames)
        if (is.null(Vcov)) Vcov <- vcov(model, dispersion = dispersion)
        if (is.vector(estimate) && all(estimate %in% seq(coefNames))) {
            if (!length(estimate))
                stop("no coefficients specified by 'estimate' argument")
            comb <- naToZero(coefs[estimate])
            var <- Vcov[estimate, estimate]
            coefMatrix <- matrix(0, l, length(comb))
            coefMatrix[cbind(estimate, seq(length(comb)))] <- 1
            colnames(coefMatrix) <- names(comb)
        }
        else {
            coefMatrix <- as.matrix(estimate)
            if (!is.numeric(coefMatrix))
                stop("'estimate' should specify parameters using ",
                     "\"pick\" or a vector of \n names/indices; ",
                     "or specify linear combinations using ",
                     "a numeric vector/matrix.")
            if (eliminate && nrow(coefMatrix) == l - eliminate)
                coefMatrix <- cbind(matrix(0, eliminate, ncol(coefMatrix)),
                                    coefMatrix)
            if (nrow(coefMatrix) != l)
                stop("NROW(estimate) should equal length(coef(model)) or \n",
                     "length(coef(model)) - model$eliminate")
            comb <- drop(crossprod(coefMatrix, naToZero(coefs)))
            var <- crossprod(coefMatrix, crossprod(Vcov, coefMatrix))
        }
    }
    estimable <- rep(TRUE, ncol(coefMatrix))
    if (checkEstimability) {
        estimable <- checkEstimable(model, coefMatrix, ...)
        if (any(!na.omit(estimable)))
            cat("Std. Error is NA where estimate is fixed or unidentified\n")
    }
    if (is.matrix(var))
        sterr <- sqrt(diag(var))
    else
        sterr <- sqrt(var)
    is.na(sterr[estimable %in% c(FALSE, NA)]) <- TRUE
    result <- data.frame(comb, sterr)
    rowNames <- colnames(coefMatrix)
    if (is.null(rowNames))
        rowNames <- paste("Combination", ncol(coefMatrix))
    dimnames(result) <- list(rowNames, c("Estimate", "Std. Error"))
    result
}
