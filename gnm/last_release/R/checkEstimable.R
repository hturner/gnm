checkEstimable <- function(model,
                            combMatrix = diag(length(coef(model))),
                            tolerance = NULL)
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
    eliminate <- model$eliminate
    if (!is.null(eliminate)) {
        ## sweeps needed to get the rank right
        subtracted <- rowsum(X, eliminate)/tabulate(eliminate)
        if (attr(terms(model), "intercept") == 1) subtracted[,1] <- 0
        X <- X - subtracted[eliminate,]
    }
    rankX <- rankMatrix(X)
    check.1 <- function(comb){
        Xc <- rbind(X, comb)
        rankXc <- rankMatrix(Xc, tol = tolerance)
        return(rankXc == rankX)
    }
    result[!resultNA] <- apply(combMatrix[, !resultNA, drop = FALSE],
                               2, check.1)
    names(result) <- colnames(combMatrix)
    return(result)
}
