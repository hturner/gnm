checkEstimable <- function(model, cmatrix, tolerance = 1e-14, verbose = TRUE){
    if (verbose) cat("checking estimability...\n")
    if (!inherits(model, "gnm")) stop("model not of class gnm")
        coefs <- coef(model)
    l <- length(coefs)
    cmatrix <- as.matrix(cmatrix)
    if (nrow(cmatrix) != l) stop(
          "cmatrix does not match coef(model)")
    Xt <- t(model.matrix(model))
    test1 <- function(cvec){
        result <- logical()
        if (sd(cvec) < tolerance) {
            is.na(result) <- TRUE
            return(result)
        }
        cvec <- (cvec - mean(cvec))/sd(cvec)
        temp <- lm(cvec ~ -1 + Xt)
        rss <- sum(residuals(temp) ^ 2)
        result <- rss < tolerance
        return(result)
    }
    apply(cmatrix, 2, test1)
}
