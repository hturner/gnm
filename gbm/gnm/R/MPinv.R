"MPinv" <-
    function (Xmat, tol = 100*.Machine$double.eps) 
{ ## Moore-Penrose pseudoinverse of a real-valued matrix.
  ## Patterned after ginv() from the MASS package of W N Venables
  ## and B D Ripley.
  ## This version retains any row and column names, and has a different
  ## default tolerance for detection of zero eigenvalues.  It also adds
  ## the computed rank as an attribute of the result.
    if (!is.matrix(Xmat)) stop("Xmat is not a matrix")
    Rownames <- rownames(Xmat)
    Colnames <- colnames(Xmat)
    Xsvd <- svd(Xmat)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    result <- {
        if (all(Positive)) 
            Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
        else if (!any(Positive)) 
            array(0, dim(Xmat)[2:1])
        else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                            t(Xsvd$u[, Positive, drop = FALSE]))
    }
    attr(result, "rank") <- sum(Positive)
    if (!is.null(Rownames)) colnames(result) <- Rownames
    if (!is.null(Colnames)) rownames(result) <- Colnames
    return(result)
}
