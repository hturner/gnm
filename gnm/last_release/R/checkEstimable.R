checkEstimable <- function(model, combMatrix = NULL,
                           tolerance = 1e6 * .Machine$double.eps)
{
    if (!inherits(model, "gnm")) stop("model not of class gnm")
    coefs <- coef(model)
    X <- model.matrix(model)
    l <- ncol(X)
    constrained <- seq(l) %in% model$constrain
    X <- X[, !constrained, drop = FALSE]
    if (is.null(combMatrix))
        combMatrix <- diag(seq_len(l))
    else combMatrix <- as.matrix(combMatrix)
    if (nrow(combMatrix) != l) stop(
          "dimensions of combMatrix do not match coef(model)")
    combMatrix <- scale(combMatrix[!constrained, ], center = FALSE)
    resultNA <- is.na(colSums(combMatrix))
    resids <- qr.resid(qr(crossprod(X)), combMatrix[, !resultNA, drop = FALSE])
    rms <- apply(resids, 2, sd)
    result <- rep.int(NA, length(resultNA))
    result[!resultNA] <- rms < tolerance
    names(result) <- colnames(combMatrix)
    return(result)
}
