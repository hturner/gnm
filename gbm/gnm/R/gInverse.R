MPinv <- function(mat, eliminate = numeric(0), first.col.only = FALSE,
                     non.elim.only = FALSE, tol = 100*.Machine$double.eps){
    ## Moore-Penrose pseudoinverse of a real-valued matrix.
    ## Patterned after ginv() from the MASS package of W N Venables
    ## and B D Ripley.
    ## This version retains any row and column names, and has a different
    ## default tolerance for detection of zero eigenvalues.  It also adds
    ## the computed rank as an attribute of the result.
    ##
    ## For symmetric matrices, there is the option of direct inversion
    ## of an assumed-diagonal submatrix T (notation as in Harville 1997,
    ## p121) to speed up the computation especially when T is large.
    ##
    if (length(eliminate) == 0){
        if (!is.matrix(mat)) stop("mat is not a matrix")
        Rownames <- rownames(mat)
        Colnames <- colnames(mat)
        Svd <- svd(mat)
        Positive <- Svd$d > max(tol * Svd$d[1], 0)
        result <- {
            if (all(Positive)) 
                Svd$v %*% (1/Svd$d * t(Svd$u))
            else if (!any(Positive)) 
                array(0, dim(mat)[2:1])
            else Svd$v[, Positive, drop = FALSE] %*% ((1/Svd$d[Positive]) * 
                                    t(Svd$u[, Positive, drop = FALSE]))
        }
        attr(result, "rank") <- sum(Positive)
        if (!is.null(Rownames)) colnames(result) <- Rownames
        if (!is.null(Colnames)) rownames(result) <- Colnames
        if (first.col.only) result <- result[, 1, drop = FALSE]
        return(result)
    }
    ## The rest is for the case length(eliminate) > 0
    if (nrow(mat) != ncol(mat)) stop(
            "mat must be a symmetric matrix" ## no more check than this!!
            )
    n <- nrow(mat)
    elim <- 1:n %in% eliminate
    diag.indices <- (n*(0:(n-1)) + 1:n)
    T <- mat[diag.indices[eliminate]]
    W <- mat[!elim, !elim, drop = FALSE]
    U <- mat[elim, !elim, drop = FALSE]
    Ti <- 1/T
    Ti.U <- Ti * U
    V.Ti <- t(Ti.U)
    Qmat <- W - crossprod(Ti.U, U)
    Qi <- MPinv(Qmat, tol = tol)
    rankQ <- attr(Qi, "rank")
    k <- length(T)
    result <- matrix(NA,
                     if (non.elim.only) n - k else n,
                     if (first.col.only) 1 else
                         if (non.elim.only) n - k else n)
    cols.notElim <- if (first.col.only) 1 else
                        if (non.elim.only) 1:(n - k) else
                           !elim
    rows.notElim <- if (non.elim.only) 1:(n - k) else !elim
    if (first.col.only) Qi <- Qi[, 1, drop = FALSE]
    result[rows.notElim, cols.notElim] <- Qi
    if (!non.elim.only){
                         temp <- - crossprod(Qi, V.Ti)
                         result[elim, cols.notElim] <- t(temp)
                     }
    if (!first.col.only && !non.elim.only){
        result[!elim, elim] <- temp
        temp <- crossprod(V.Ti, Qi) %*% V.Ti
        diag.indices <- k*(0:(k-1)) + 1:k
        temp[diag.indices] <- Ti + temp[diag.indices]
        result[elim, elim] <- temp
    }
    attr(result, "rank") <- rankQ + k
    theNames <- colnames(mat)
    rownames(result) <- if (non.elim.only) theNames[!elim]
                        else theNames
    colnames(result) <-
        if (first.col.only) theNames[!elim][1]
        else if (non.elim.only) theNames[!elim]
             else theNames
    result
}
