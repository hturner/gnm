checkEstimable <- function(model, coefMatrix, tolerance = 1e-8){
    if (!inherits(model, "gnm")) stop("model not of class gnm")
    coefs <- coef(model)
    l <- length(coefs)
    coefMatrix <- as.matrix(coefMatrix)
    if (nrow(coefMatrix) != l) stop(
          "coefMatrix does not match coef(model)")
    Xt <- t(model.matrix(model))
    Xt <- Xt[!is.na(coefs), ]
    coefMatrix <- coefMatrix[!is.na(coefs), ]
    coefMatrix <- scale(coefMatrix, center = FALSE)
    resultNA <- apply(coefMatrix, 2, function(col) any(is.na(col)))
    result <- logical(ncol(coefMatrix))
    is.na(result) <- resultNA
    resids <- qr.resid(qr(Xt), coefMatrix[, !resultNA, drop = FALSE])
    rss <- apply(resids, 2, var)
    result[!resultNA] <- rss < tolerance
    names(result) <- colnames(coefMatrix)
    return(result)
}
