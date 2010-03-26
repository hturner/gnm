solve1 <- function (W, Tvec = NULL, U = NULL, elim = NULL, scale = 1,
                    elimscale = 1)
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
    I1 <- numeric(nrow(W))
    I1[1] <- 1
    if (is.null(Tvec)) { ## the basic routine, no eliminated submatrix
        Qi <- solve(W, I1, tol = .Machine$double.eps)/scale
        return(-Qi[-1]/Qi[1])
    }
##  Now allow for the possibility of an eliminated submatrix
    Ti.U <- U/Tvec
    Qi <- solve(W - crossprod(U, Ti.U), I1, tol = .Machine$double.eps)
    structure(-Qi[-1]/Qi[1] * scale[1]/scale[-1],
              eliminated = (Ti.U %*% Qi)/Qi[1] * scale[1]/elimscale)
}
