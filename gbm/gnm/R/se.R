se <- function(model, cmatrix, check.identifiability = TRUE, ...){
    if (!inherits(model, "gnm")) stop(
                "model is not of class \"gnm\"")
    if (!is.numeric(cmatrix)) stop("cmatrix not numeric")
    cmatrix <- as.matrix(cmatrix)
    coefs <- coef(model)
    l <- length(coefs)
    identifiable <- rep(TRUE, ncol(cmatrix)) 
    if (nrow(cmatrix) != l) stop(
          "nrow(cmatrix) does not match length(coef(model))")
    comb <- drop(crossprod(cmatrix, coefs))
    if (check.identifiability){
        identifiable <- checkIdentifiability(model, cmatrix, ...)
    }
    var <- crossprod(cmatrix, crossprod(vcov(model), cmatrix))
    sterr <- sqrt(diag(var))
    is.na(sterr[!identifiable]) <- is.na(comb[!identifiable]) <- TRUE
    data.frame(estimate = comb, se = sterr, row.names = colnames(cmatrix))
}
