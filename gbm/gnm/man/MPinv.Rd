\name{MPinv}
\alias{MPinv}
\title{ Moore-Penrose Pseudoinverse of a Real-valued Matrix }
\description{
  Computes the Moore-Penrose generalized inverse, optionally using a
  method (appropriate for symmetric matrices only) which exploits the
  direct inversion of a nonsingular diagonal submatrix if such exists. 
}
\usage{
MPinv(mat, eliminate = numeric(0), onlyFirstCol = FALSE,
      onlyNonElim = FALSE, tolerance = 100*.Machine$double.eps,
      rank = NULL, method = "svd")
}
\arguments{
  \item{mat}{ A real matrix, which should be symmetric unless
    \code{eliminate} has zero length. }
  \item{eliminate}{ Numeric.  A vector of indices identifying
    rows/columns which form a diagonal submatrix of \code{mat}. }
  \item{onlyFirstCol}{ Logical.  If \code{TRUE}, return only the first
    non-eliminated column of the generalized inverse matrix; otherwise
    the whole matrix.}
  \item{onlyNonElim}{ Logical.  If \code{TRUE}, return only rows/columns that
    are in the complement of the set specified in \code{eliminate}.}
  \item{tolerance}{ A positive scalar which determines the tolerance for
    detecting zeroes among the singular values. }
  \item{rank}{Either \code{NULL}, in which case the rank of \code{mat} is
    determined numerically; or an integer specifying the rank of
    \code{mat} if it is known.  No check is made on the validity of any
    non-\code{NULL} value.}
  \item{method}{Character, one of \code{"svd", "chol"}}
}
\details{
  Real-valuedness is not checked, neither is symmetry, nor the
  diagonality of the submatrix when \code{eliminate} is used.

  When \code{eliminate} is used, all of the eliminated submatrix's
  diagonal entries must be non-zero.  

  The purpose of the \code{eliminate} argument is a substantial
  reduction of the computational burden when the number of `eliminated'
  rows/columns is large.

  The specification \code{method = "chol"} is valid only for
  symmetric matrices.  No check for symmetry is made.
}
\value{
  A matrix, with an additional attribute named \code{"rank"} containing
  the numerically determined rank of the matrix.
}
\references{Harville, D. A. (1997).  \emph{Matrix Algebra from a
    Statistician's Perspective}.  New York: Springer.  }
\author{ David Firth }
\seealso{\code{\link{ginv}}}
\examples{
A <- matrix(c(1, 1, 0,
              1, 1, 0,
              2, 3, 4), 3, 3)
B <- MPinv(A)
A \%*\% B \%*\% A - A  # essentially zero
B \%*\% A \%*\% B - B  # essentially zero
attr(B, "rank")    # here 2

## demonstration that "svd" and "chol" deliver essentially the same
## results for symmetric matrices:
A <- crossprod(A)
MPinv(A) - MPinv(A, method = "chol") ##  (essentially zero) 

## now a symmetric example with diagonal submatrix to `eliminate'
A <- matrix(c(1, 0, 2,
              0, 2, 3,
              2, 3, 4), 3, 3)
B <- MPinv(A, eliminate = 1:2)
A \%*\% B \%*\% A - A  # essentially zero
B \%*\% A \%*\% B - B  # essentially zero
attr(B, "rank")        # here 3  

\dontrun{
## demo that eliminate can give substantial speed gains
A <- diag(rnorm(100))
A <- cbind(A, matrix(rnorm(200), 100, 2))
A <- rbind(A, cbind(t(A[, 101:102]), matrix(c(1, 2, 2, 1), 2, 2)))
system.time(for (i in 1:1000) B <-  MPinv(A))
##  [1] 30.85  0.06 30.91  0.00  0.00
system.time(for (i in 1:1000) B <-  MPinv(A, eliminate = 1:100))
##  [1] 3.10 0.00 3.11 0.00 0.00
}
}
\keyword{ array }

