se <- function(model, estimate = "all", checkEstimability = TRUE, ...){
    if (!inherits(model, "gnm")) stop("model is not of class \"gnm\"")
    coefs <- coef(model)
    l <- length(coefs)
    eliminate <- model$eliminate
    if (eliminate && eliminate == l)
        stop("No non-eliminated coefficients")
    coefNames <- names(coefs)[(eliminate + 1):l]
    if (identical(estimate, "all")) {
        if (eliminate == 0)
            return(data.frame(coef(summary(model)))[, 1:2])
        return(data.frame(coef(summary(model)))[(eliminate + 1):l, 1:2])
    }
    else {
        if (identical(estimate, "pick")) {
            estimate <-
                unlist(relimp:::pickFrom(coefNames,
                                         setlabels = "Selected coefficients",
                                         title =
                                         paste("Estimate standard errors",
                                         "for one or more gnm coefficients"),
                                         items.label = "Model coefficients:",
                                         edit.setlabels = FALSE))
            if(!length(nchar(estimate)))
                stop("no parameters were selected")
        }
        if (is.character(estimate))
            estimate <- match(estimate, coefNames)
        if (is.vector(estimate) && all(estimate %in% seq(l - eliminate))) {
            if (!length(estimate))
                stop("no coefficients specified by 'estimate' argument")
            estimate <- estimate + eliminate
            comb <- naToZero(coefs[estimate])
            var <- vcov(model)[estimate, estimate]
            coefMatrix <- matrix(0, l, length(comb))
            coefMatrix[cbind(estimate, seq(comb))] <- 1
            colnames(coefMatrix) <- names(comb)
        }
        else {
            coefMatrix <- as.matrix(estimate)
            if (!is.numeric(coefMatrix))
                stop("'estimate' must specify parameters using ",
                     "\"all\", \"pick\", or a vector of \n names/indices; ",
                     "or specify linear combinations using ",
                     "a numeric vector/matrix.")
            if (eliminate && nrow(coefMatrix) == l - eliminate)
                coefMatrix <- cbind(matrix(0, eliminate, ncol(coefMatrix)),
                                    coefMatrix)                
            if (nrow(coefMatrix) != l)
                stop("NROW(estimate) should equal length(coef(model)) or \n",
                     "length(coef(model)) - model$eliminate")
            comb <- drop(crossprod(coefMatrix, naToZero(coefs)))
            var <- crossprod(coefMatrix, crossprod(vcov(model), coefMatrix))
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
