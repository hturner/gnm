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
    n <- nrow(W)
    if (is.null(Tvec)) { ## the basic routine, no eliminated submatrix
        I1 <- numeric(n)
        I1[1] <- 1
        return(drop(solve(W, I1)))
    }
##  Now allow for the possibility of an eliminated submatrix
    Ti <- sqrt(1/Tvec)
    Ti.U <- Ti * U
    Qi <- solve1(W - crossprod(Ti.U))
    result <- numeric(n + length(Tvec) - 1)
    result[elim] <- -tcrossprod(Qi, Ti * Ti.U)
    result[-elim] <- Qi[-1]
    -result/Qi[1]
}
