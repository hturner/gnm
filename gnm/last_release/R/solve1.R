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
        return(drop(solve(W, I1, tol = .Machine$double.eps)))
    }
##  Now allow for the possibility of an eliminated submatrix
    Ti.U <- U/Tvec
    Qi <- solve1(W - crossprod(U, Ti.U))
    result <- list(coefficients = -Qi[-1]/Qi[1],
                   elim.coefs = (Ti.U %*% Qi)/Qi[1])
}
