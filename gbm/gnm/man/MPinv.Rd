\name{MPinv}
\alias{MPinv}
\title{ Moore-Penrose inverse of a real-valued matrix }
\description{
  Computes the Moore-Penrose generalized inverse (or pseudoinverse)
  of a real-valued matrix, using the singular value decomposition.
}
\usage{
MPinv(Xmat, tol = 100 * .Machine$double.eps)
}
\arguments{
  \item{Xmat}{ A real-valued matrix }
  \item{tol}{ A positive scalar which determines the tolerance for
  detecting zeroes among the singular values }
}
\value{
  A matrix.
}
\author{ David Firth }
\note{ Similar to \code{\link{MASS:ginv}}, but without the facility
to handle complex-valued matrices, and with a different default tolerance.
\code{MPinv} also preserves row and column names where available, and 
returns the rank of the matrix as an attribute of the result. }
\examples{
A <- matrix(c(
 1, 1, 0,
 1, 1, 0,
 2, 3, 4), 3, 3)
B <- MPinv(A)
B %*% A %*% B - B  ## essentially zero
attr(B, "rank")    ## here 2
}
\keyword{ array }
