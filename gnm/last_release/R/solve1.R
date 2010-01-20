solve1 <- function (W, Tvec = NULL, U = NULL, elim = NULL)
##  Inverse of partitioned matrix, essentially
##    T  U
##    U' W
##  where D is diag(Tvec).
##
##  ("essentially", because T does not necessarily occupy the *first*
##   k rows and columns; the `elim' argument specifies which rows/columns of
##   the matrix are actually occupied by T.)
##
##  Result is only the first row/column of the inverse
{
    if (is.null(Tvec)) { ## the basic routine, no eliminated submatrix
        I1 <- numeric(nrow(W))
        I1[1] <- 1
        return(drop(solve(W, I1)))
    }
##  Now allow for the possibility of an eliminated submatrix
    n <- ncol(W)
    Ti <- sqrt(1/Tvec)
    k <- length(Tvec)
    Ti.U <- Ti * U
    Qmat <- W - crossprod(Ti.U)
    Qi <- solve1(Qmat)
    result <- numeric(n + k)
    result[!elim] <- Qi
    result[elim] <- -tcrossprod(Qi, Ti * Ti.U)
    result
}
