qrSolve <- function(A, b, rank = NULL, ...) {
    ## Function to solve Ax = b, giving the minimum norm solution
    ## (which corresponds to use of the Moore-Penrose inverse) in
    ## the non-full-rank case
    theRank <- rank
    mA <- nrow(A)
    nA <- ncol(A)
    bMissing <- FALSE
    if (missing(b)){
        if (nA != mA) stop("only square matrices can be inverted")
        b <- diag(nA)
        bMissing <- TRUE
    }
    widthA <- min(mA, nA)
    qrA <- qr(A, ...)
    rA <- if (is.null(theRank)) qrA$rank else theRank
    ## Full-rank case:
    if (rA == widthA) return(solve(qrA, b))
    ## Non-full-rank case: using the complete orthogonal decomposition
    ## (see for example lecture notes by James P Reilly, at
    ##    http://www.ece.mcmaster.ca/~reilly/ee731/ch10.ps
    ## )
    Rtrans <- t(qr.R(qrA)[1:rA, ]) ## foot of Reilly's p7
    reversed <- rA:1
    qr2 <- qr(Rtrans[, reversed], ...)
    Ttrans <- qr.R(qr2)[reversed, reversed] ## so it's lower triangular
    Ztrans <- t(qr.Q(qr2)[, reversed])
    cmat <- crossprod(qr.Q(qrA)[, 1:rA], b)
    w <- forwardsolve(Ttrans, cmat, transpose = TRUE)
    result <- crossprod(Ztrans, w)[order(qrA$pivot), ]
    attr(result, "rank") <- rA
    if (is.matrix(result)) {
        rownames(result) <- colnames(A)
        colnames(result) <- if (bMissing) rownames(A)
        else if (is.matrix(b)) colnames(b)
        else names(b)
    } else names(result) <- colnames(A)
    result
}
