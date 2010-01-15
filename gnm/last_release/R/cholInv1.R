cholInv1 <- function (W, Tvec = NULL, U = NULL, elim = NULL,
                      result.as.vector = TRUE)
##  Inverse of partitioned matrix, essentially
##    T  U
##    U' W
##  where D is diag(Tvec).
##
##  ("essentially", because T does not necessarily occupy the *first*
##   k rows and columns; the `elim' argument specifies which rows/columns of
##   the matrix are actually occupied by T.)
##
##  Result is only the first row of the inverse, either as a vector or as a
##  one-row matrix.
{
    if (is.null(Tvec)) { ## the basic routine, no eliminated submatrix
        result <- chol2inv(chol(W))
        result <- result[1, , drop = result.as.vector]
        return(result)
    }
##  Now allow for the possibility of an eliminated submatrix
    n <- ncol(W)
    Ti <- sqrt(1/Tvec)
    k <- length(Tvec)
    elim <- {if (is.null(elim)) c(rep(TRUE, k), rep(FALSE, n))
             else seq(n + k) %in% elim}
    nonElim <- which(!elim)
    Ti.U <- Ti * U
    Qmat <- W - crossprod(Ti.U)
    Qi <- as.matrix(cholInv1(Qmat, result.as.vector = FALSE)) #else S4 Matrix class
    result <- numeric(n + k)
    result[nonElim] <- Qi
    result[elim] <- -tcrossprod(Ti * Ti.U, Qi)
    if (result.as.vector) result else matrix(result, 1, n + k)
}
