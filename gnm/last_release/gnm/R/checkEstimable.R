checkEstimable <- function(model, combMatrix = diag(seq(along = coef(model))),
                           tolerance = 1e6 * .Machine$double.eps)
{
    if (!inherits(model, "gnm")) stop("model not of class gnm")
    coefs <- coef(model)
    l <- length(coefs)
    combMatrix <- as.matrix(combMatrix)
    if (nrow(combMatrix) != l) stop(
          "dimensions of combMatrix do not match coef(model)")
    X <- model.matrix(model)[, !is.na(coefs), drop = FALSE]
    combMatrix <- scale(combMatrix[!is.na(coefs), ], center = FALSE)
    resultNA <- apply(combMatrix, 2, function(col) any(is.na(col)))
    result <- logical(ncol(combMatrix))
    is.na(result) <- resultNA
    resids <- qr.resid(qr(crossprod(X)), combMatrix[, !resultNA, drop = FALSE])
    rms <- apply(resids, 2, sd)
    result[!resultNA] <- rms < tolerance
    names(result) <- colnames(combMatrix)
    return(result)
}
