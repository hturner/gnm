se <- function(model, estimate = "all", checkEstimability = TRUE, ...){
    if (!inherits(model, "gnm")) stop("model is not of class \"gnm\"")
    coefs <- coef(model)
    coefs[is.na(coefs)] <- 0
    coefNames <- names(coefs)
    l <- length(coefs)
    if (identical(estimate, "all"))
        coefMatrix <- diag(l)
    else {
        if (identical(estimate, "pick")) {
            if (!require(tcltk))
                stop("estimate = \"pick\", and tcltk not installed")
            if (!require(relimp))
                stop("the relimp package from CRAN needs to be installed")
            estimate <-
                unlist(pickFrom(coefNames[(model$eliminate + 1):l],
                                setlabels = "Selected coefficients",
                                title = paste("Estimate standard errors",
                                "for one or more gnm coefficients"),
                                items.label = "Model coefficients:",
                                edit.setlabels = FALSE))
            if(!length(nchar(estimate)))
                stop("no parameters were selected")
        }
        if (is.character(estimate))
            estimate <- match(estimate, coefNames)
        if (is.vector(estimate)) {
            if (!length(estimate))
                stop("no coefficients specified by 'estimate' argument")
            estimate <- cbind(estimate, seq(along = estimate))
            coefMatrix <- matrix(0, l, nrow(estimate))
            coefMatrix[estimate] <- 1
            colnames(coefMatrix) <- coefNames[estimate[,1]]
        }
        else {
            coefMatrix <- as.matrix(estimate)
            if (!is.numeric(coefMatrix))
                stop("'estimate' argument must be one of \"all\", \"pick\", a ",
                     "character vector of names \n or a matrix of linear ",
                     "combinations")
            if (nrow(coefMatrix) != l)
                stop("nrow(estimate) does not match length(coef(model))")
        }
    }
    estimable <- rep(TRUE, ncol(coefMatrix))
    comb <- drop(crossprod(coefMatrix, coefs))
    if (checkEstimability){
        estimable <- checkEstimable(model, coefMatrix, ...)
        if (any(is.na(estimable)))
            warning("Not all of the desired estimates are identifiable")
    }
    var <- crossprod(coefMatrix, crossprod(vcov(model), coefMatrix))
    sterr <- sqrt(diag(var))
    is.na(sterr[!estimable]) <- is.na(comb[!estimable]) <- TRUE
    data.frame(estimate = comb, SE = sterr, row.names = colnames(coefMatrix))
}
