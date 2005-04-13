checkEstimable <- function(model, cmatrix, tolerance = 1e-8){
    if (!inherits(model, "gnm")) stop("model not of class gnm")
        coefs <- coef(model)
    l <- length(coefs)
    cmatrix <- as.matrix(cmatrix)
    if (nrow(cmatrix) != l) stop(
          "cmatrix does not match coef(model)")
    if (model$eliminate > 0) {
        extra.rows <- matrix(0, model$eliminate, ncol(cmatrix))
        cmatrix <- rbind(extra.rows, cmatrix)
    }
    Xt <- t(model.matrix(model))
    cmatrix <- scale(cmatrix)
    resultNA <- apply(cmatrix, 2, function(col) any(is.na(col)))
    result <- logical(ncol(cmatrix))
    is.na(result) <- resultNA
    resids <- qr.resid(qr(Xt), cmatrix[, !resultNA, drop = FALSE]) 
    rss <- apply(resids, 2, var)
    result[!resultNA] <- rss < tolerance
    names(result) <- colnames(cmatrix)
    return(result)
}
