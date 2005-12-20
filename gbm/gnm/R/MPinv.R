MPinv <- function(mat, eliminate = numeric(0), onlyFirstCol = FALSE,
                     onlyNonElim = FALSE, tolerance = 100*.Machine$double.eps,
                     rank = NULL){
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
        if (is.null(rank)) {
            Positive <- Svd$d > max(tolerance * Svd$d[1], 0)
        } else Positive <- c(rep(TRUE, rank),
                             rep(FALSE, length(Svd$d) - rank))
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
        if (onlyFirstCol) result <- result[, 1, drop = FALSE]
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
    if (any(T == 0)) stop(
         "an eliminated submatrix must have all its diagonal entries non-zero."
         )
    Ti <- 1/T
    ## Special case is when the whole matrix is "eliminated"
    if (length(eliminate) == nrow(mat)) {
        if (onlyFirstCol) {
            result <- matrix(c(Ti[1], rep(0, n-1)), n, 1)
        } else if (onlyNonElim) {
            result <- matrix(numeric(0), 0, 0)
        } else result <- diag(Ti)
        attr(result, "rank") <- n
    } else {
        W <- mat[!elim, !elim, drop = FALSE]
        U <- mat[elim, !elim, drop = FALSE]
        Ti.U <- Ti * U
        V.Ti <- t(Ti.U)
        Qmat <- W - crossprod(Ti.U, U)
        Qi <- MPinv(Qmat, tolerance = tolerance)
        rankQ <- attr(Qi, "rank")
        k <- length(T)
        result <- matrix(NA,
                         if (onlyNonElim) n - k else n,
                         if (onlyFirstCol) 1 else
                         if (onlyNonElim) n - k else n)
        cols.notElim <- if (onlyFirstCol) 1 else {
            if (onlyNonElim) 1:(n - k) else !elim}
        rows.notElim <- if (onlyNonElim) 1:(n - k) else !elim
        if (onlyFirstCol) Qi <- Qi[, 1, drop = FALSE]
        result[rows.notElim, cols.notElim] <- Qi
        if (!onlyNonElim){
            temp <- - crossprod(Qi, V.Ti)
            result[elim, cols.notElim] <- t(temp)
        }
        if (!onlyFirstCol && !onlyNonElim){
            result[!elim, elim] <- temp
            temp <- crossprod(V.Ti, Qi) %*% V.Ti
            diag.indices <- k*(0:(k-1)) + 1:k
            temp[diag.indices] <- Ti + temp[diag.indices]
            result[elim, elim] <- temp
        }
        attr(result, "rank") <- rankQ + k
    }
    theNames <- colnames(mat)
    if (length(result) > 0.5) { ## matrix is not 0 by 0
        rownames(result) <- if (onlyNonElim) theNames[!elim]
                            else theNames
        colnames(result) <-
            if (onlyFirstCol) theNames[!elim][1]
            else if (onlyNonElim) theNames[!elim]
            else theNames
    }
    return(result)
}
