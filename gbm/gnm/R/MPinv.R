MPinv <- function (mat,
                   eliminate = numeric(0),
                   onlyFirstCol = FALSE,
                   onlyNonElim = FALSE,
                   tolerance = 100 * .Machine$double.eps,
                   rank = NULL,
                   method = "svd")
{
    theRank <- rank
    m <- nrow(mat)
    n <- ncol(mat)
    if (length(eliminate) == 0) { ## the basic routine, no eliminated submatrix
        if (!is.matrix(mat))
            stop("mat is not a matrix")
        Rownames <- rownames(mat)
        Colnames <- colnames(mat)
        if (method == "svd") {
            Svd <- svd(mat)
            Positive <- rep(FALSE, length(Svd$d))
            if (is.null(theRank)) {
                Positive <- Svd$d > max(tolerance * Svd$d[1], 0)
            } else Positive[1:theRank] <- TRUE
            result <- {
                if (all(Positive))
                    Svd$v %*% (1/Svd$d * t(Svd$u))
                else if (!any(Positive))
                    array(0, dim(mat)[2:1])
                else Svd$v[, Positive, drop = FALSE] %*% ((1/Svd$d[Positive]) *
                                       t(Svd$u[, Positive, drop = FALSE]))
            }
            attr(result, "rank") <- sum(Positive)
        }
        if (method == "chol") {
            ## Generalized inverse of a symmetric matrix using a
            ## streamlined version of the "fast" method of
            ## Courrieu, P. (2005).  Fast computation of Moore-Penrose
            ## inverse matrices. Neural Information Processing 8, 25-29.
            ##
            ## No test for symmetry performed here!
            if (!(m == n)) stop("the matrix is not symmetric")
            S <- chol(mat, pivot = TRUE) ## non-full-rank case
            if (is.null(theRank)) {
                theRank <- qr(S)$rank ## fails only on the bwt.po example
               # theRank <- attr(S, "rank") ## seems less reliable in general
            }
            pivot <- attr(S, "pivot")
            oPivot <- order(pivot)
            Lt <- S[oPivot[oPivot %in% 1:theRank], oPivot]
            L <- t(Lt)
            LLinv <- chol2inv(chol(crossprod(L)))
            result <- crossprod(Lt, crossprod(LLinv)) %*% Lt
            attr(result, "rank") <- theRank
        }
        if (!is.null(Rownames))
            colnames(result) <- Rownames
        if (!is.null(Colnames))
            rownames(result) <- Colnames
        if (onlyFirstCol)
            result <- result[, 1, drop = FALSE]
        return(result)
    }
##  Now allow for the possibility of an eliminated submatrix
    if (m != n)
        stop("mat must be a symmetric matrix")
    n <- nrow(mat)
    elim <- 1:n %in% eliminate
    diag.indices <- (n * (0:(n - 1)) + 1:n)
    T <- mat[diag.indices[eliminate]]
    if (any(T == 0))
      stop("an eliminated submatrix must have all diagonal entries non-zero.")
    if (all(elim)) {
    ## Are *all* rows/cols eliminated (ie matrix is diagonal)?
        if (onlyNonElim || onlyFirstCol) {
            stop("there are no non-eliminated rows/columns")
        }
        result <- diag(1/T)
        attr(result, "rank") <- m
    }
    else {
        W <- mat[!elim, !elim, drop = FALSE]
        U <- mat[elim, !elim, drop = FALSE]
        Ti <- 1/T
        k <- length(T)
        Ti.U <- Ti * U
        V.Ti <- t(Ti.U)
        Qmat <- W - crossprod(Ti.U, U)
        rankQ <- if (is.null(theRank)) NULL else theRank - length(T)
        Qi <- MPinv(Qmat, tolerance = tolerance, rank = rankQ,
                    method = method)
        rankQ <- attr(Qi, "rank")
        result <- matrix(NA,
                         if (onlyNonElim) n - k else n,
                         if (onlyFirstCol) 1 else {
                             if (onlyNonElim) n - k else n
                             }
                         )
        cols.notElim <- if (onlyFirstCol) 1 else {
            if (onlyNonElim) 1:(n - k) else !elim
        }
        rows.notElim <- if (onlyNonElim) 1:(n - k) else !elim
        if (onlyFirstCol)
            Qi <- Qi[, 1, drop = FALSE]
        result[rows.notElim, cols.notElim] <- Qi
        if (!onlyNonElim) {
            temp <- -crossprod(Qi, V.Ti)
            result[elim, cols.notElim] <- t(temp)
        }
        if (!onlyFirstCol && !onlyNonElim) {
            result[!elim, elim] <- temp
            temp <- crossprod(V.Ti, Qi) %*% V.Ti
            diag.indices <- k * (0:(k - 1)) + 1:k
            temp[diag.indices] <- Ti + temp[diag.indices]
            result[elim, elim] <- temp
        }
        attr(result, "rank") <- rankQ + k
    }
    theNames <- colnames(mat)
    rownames(result) <- if (onlyNonElim) theNames[!elim] else theNames
    colnames(result) <- if (onlyFirstCol) theNames[!elim][1] else {
        if (onlyNonElim) theNames[!elim] else theNames
    }
    result
}
