MPinv <- function (mat,
                   tolerance = 100 * .Machine$double.eps,
                   rank = NULL,
                   method = "svd")
{
    theRank <- rank
    if (!is.matrix(mat)) stop("mat must be a matrix")
    m <- nrow(mat)
    n <- ncol(mat)
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
            else Svd$v[, Positive, drop = FALSE] %*%
                ((1/Svd$d[Positive]) *
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
        S <- suppressWarnings(chol(mat, pivot = TRUE))
        ## (non-full-rank case)
        if (is.null(theRank)) {
            theRank <- qr(S)$rank ## fails only on the bwt.po example
            ## theRank <- attr(S, "rank") ## seems less reliable in general
        }
        pivot <- attr(S, "pivot")
        oPivot <- order(pivot)
        Lt <- S[oPivot[oPivot %in% 1:theRank], oPivot, drop = FALSE]
        LLinv <- chol2inv(chol(tcrossprod(Lt)))
        result <- crossprod(Lt, crossprod(LLinv)) %*% Lt
        attr(result, "rank") <- theRank
    }
    if (!is.null(Rownames))
        colnames(result) <- Rownames
    if (!is.null(Colnames))
        rownames(result) <- Colnames
    return(result)
}
