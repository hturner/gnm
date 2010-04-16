\name{MPinv}
\alias{MPinv}
\title{ Moore-Penrose Pseudoinverse of a Real-valued Matrix }
\description{
  Computes the Moore-Penrose generalized inverse. 
}
\usage{
MPinv(mat, tolerance = 100*.Machine$double.eps,
      rank = NULL, method = "svd")
}
\arguments{
  \item{mat}{ a real matrix.}
  \item{tolerance}{ A positive scalar which determines the tolerance for
    detecting zeroes among the singular values. }
  \item{rank}{Either \code{NULL}, in which case the rank of \code{mat} is
    determined numerically; or an integer specifying the rank of
    \code{mat} if it is known.  No check is made on the validity of any
    non-\code{NULL} value.}
  \item{method}{Character, one of \code{"svd", "chol"}. The
    specification \code{method = "chol"} is valid only for 
    symmetric matrices.  }
}
\details{
  Real-valuedness is not checked, neither is symmetry when \code{method
    = "chol"}. 
}
\value{
  A matrix, with an additional attribute named \code{"rank"} containing
  the numerically determined rank of the matrix.
}
\references{
  Harville, D. A. (1997).  \emph{Matrix Algebra from a
    Statistician's Perspective}.  New York: Springer.

  Courrieu, P. (2005).  Fast computation of Moore-Penrose
  inverse matrices. \emph{Neural Information Processing} \bold{8},
  25--29
}
\author{ David Firth }
\seealso{\code{\link[MASS]{ginv}}}
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
}
\keyword{ array }

