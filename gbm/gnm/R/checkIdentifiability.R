checkIdentifiability <-
    function(model, cmatrix, tolerance = 1e-4, verbose = TRUE){
    if (verbose) cat("checking identifiability...\n")
    if (!inherits(model, "gnm")) stop("model not of class gnm")
    coefs <- coef(model)
        l <- length(coefs)
    cmatrix <- as.matrix(cmatrix)
    if (nrow(cmatrix) != l) stop(
          "cmatrix does not match coef(model)")
    estimates <- drop(crossprod(cmatrix, coefs))
    temp <- update(model, trace = FALSE, start = NULL)
    estimates2 <- drop(crossprod(cmatrix, coef(temp)))
    result <- abs(estimates2 - estimates) < tolerance
    names(result) <- colnames(cmatrix)
    return(result)
}

    
