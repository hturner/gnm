qrSolve <- function(A, b, theRank = NULL) {
    ## Function to solve Ax = b, giving the minimum norm solution
    ## (which corresponds to use of the Moore-Penrose inverse) in
    ## the non-full-rank case
    mA <- nrow(A)
    nA <- ncol(A)
    widthA <- min(mA, nA)
    qrA <- qr(A)
    rA <- if (is.null(theRank)) qrA$rank else theRank
    ## Full-rank case:
    rA <- qrA$rank
    if (rA == widthA) return(solve(qrA, b))
    ## Non-full-rank case: using the complete orthogonal decomposition
    ## (see for example lecture notes by James P Reilly, at
    ##    http://www.ece.mcmaster.ca/~reilly/ee731/ch10.ps
    ## )
    ##
    ## But this is acutally slower (because qr is slower than chol)
    ## for symmetric A, than is MPinv(..., method = "chol")
    Rtrans <- t(qr.R(qrA)[1:rA, ]) ## foot of Reilly's p7
    reversed <- rA:1
    qr2 <- qr(Rtrans[, reversed])
    Ttrans <- qr.R(qr2)[reversed, reversed] ## so it's lower triangular
  #  Ttrans <- qr.R(qr2)
    Ztrans <- t(qr.Q(qr2)[, reversed])
    c <- crossprod(qr.Q(qrA)[, 1:rA], b)
    w <- forwardsolve(Ttrans, c, transpose = TRUE)
    y <- rep(0, widthA - rA)
    crossprod(Ztrans, w)[order(qrA$pivot)]
}
