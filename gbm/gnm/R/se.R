se <- function(model, coefMatrix, checkEstimability = TRUE, ...){
    if (!inherits(model, "gnm")) stop("model is not of class \"gnm\"")
    if (!is.numeric(coefMatrix)) stop("coefMatrix not numeric")
    coefMatrix <- as.matrix(coefMatrix)
    coefs <- coef(model)
    l <- length(coefs)
    estimable <- rep(TRUE, ncol(coefMatrix))
    if (nrow(coefMatrix) != l) stop(
          "nrow(coefMatrix) does not match length(coef(model))")
    comb <- drop(crossprod(coefMatrix, coefs))
    if (checkEstimability){
        estimable <- checkEstimable(model, coefMatrix, ...)
    }
    var <- crossprod(coefMatrix, crossprod(vcov(model), coefMatrix))
    sterr <- sqrt(diag(var))
    is.na(sterr[!estimable]) <- is.na(comb[!estimable]) <- TRUE
    data.frame(estimate = comb, SE = sterr, row.names = colnames(coefMatrix))
}
