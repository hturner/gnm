se <- function(model, cmatrix, check.estimability = TRUE, ...){
    if (!inherits(model, "gnm")) stop(
                "model is not of class \"gnm\"")
    if (!is.numeric(cmatrix)) stop("cmatrix not numeric")
    cmatrix <- as.matrix(cmatrix)
    coefs <- coef(model)
    l <- length(coefs)
    estimable <- rep(TRUE, ncol(cmatrix)) 
    if (nrow(cmatrix) != l) stop(
          "nrow(cmatrix) does not match length(coef(model))")
    comb <- drop(crossprod(cmatrix, coefs))
    if (check.estimability){
        estimable <- checkEstimable(model, cmatrix, ...)
    }
    var <- crossprod(cmatrix, crossprod(vcov(model), cmatrix))
    sterr <- sqrt(diag(var))
    is.na(sterr[!estimable]) <- is.na(comb[!estimable]) <- TRUE
    data.frame(estimate = comb, se = sterr, row.names = colnames(cmatrix))
}
