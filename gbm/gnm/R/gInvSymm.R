gInvSymm <- function(mat, eliminate = numeric(0), first.col.only = FALSE,
                     non.elim.only = FALSE, ...){
    ## Notation as in Harville, "Matrix Algebra from a Statistician's
    ## Perspective", page 121
    ##
    ## G-inverse with direct inversion of an assumed-diagonal submatrix T
    if (length(eliminate) == 0)
        return(if (first.col.only) MPinv(mat, ...)[, 1, drop = FALSE]
               else MPinv(mat, ...)
               )
    if (nrow(mat) != ncol(mat)) stop(
            "mat must be a symmetric matrix"
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
    Qi <- MPinv(Qmat)
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
