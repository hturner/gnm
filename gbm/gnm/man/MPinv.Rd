\name{MPinv}
\alias{MPinv}
\title{ Moore-Penrose pseudoinverse of a real-valued matrix }
\description{
  Computes the Moore-penrose generalized inverse, optionally using a
  method (appropriate for symmetric matrices only) which exploits the
  direct inversion of a diagonal submatrix if such exists. 
}
\usage{
MPinv(mat, eliminate = numeric(0), first.col.only = FALSE,
non.elim.only = FALSE, tol = 100*.Machine$double.eps)
}
\arguments{
  \item{mat}{ A real matrix, which should be symmetric unless
    \code{eliminate} has zero length. }
  \item{eliminate}{ Numeric.  A vector of indices identifying
    rows/columns which form a diagonal submatrix of \code{mat}. }
  \item{first.col.only}{ Logical.  If \code{TRUE}, return only the first
    column of the generalized inverse matrix; otherwise the whole matrix.}
  \item{non.elim.only}{ Logical.  If TRUE, return only rows/columns that
    are in the complement of the set specified in \code{eliminate}.}
  \item{tol}{ A positive scalar which determines the tolerance for
  detecting zeroes among the singular values }
}
\details{
  Real-valuedness is not checked, and neither is symmetry or the
  diagonality of the submatrix when \code{eliminate} is used.

  The purpose of the \code{eliminate} argument is a substantial
  reduction of the computational burden when the number of `eliminated'
  rows/columns is large.  
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
##  [1] 26.83 10.24 39.76  0.00  0.00
system.time(for (i in 1:1000) B <-  MPinv(A, eliminate = 1:100))
##  [1] 3.49 1.37 5.04 0.00 0.00
}
}
\keyword{ array }

