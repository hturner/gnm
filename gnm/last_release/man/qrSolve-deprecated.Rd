\name{qrSolve-deprecated}
\alias{qrSolve-deprecated}
\title{ Minimum-length solution of a linear system }
\description{
  This function solves the linear system Ax = b, producing
  the minimum-length solution if \code{A} is not of full rank.
  If b is not specified, the Moore-Penrose generalized inverse
  (pseudo-inverse) is returned.
  The method involves two applications of QR decomposition.
}
\usage{
qrSolve(A, b, rank = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ A numeric matrix (m by n, say) }
  \item{b}{ A numeric vector of length m, or a matrix with m rows }
  \item{rank}{ Numeric: the algebraic rank of \code{A}, if known. }
  \item{\dots}{ Other possible arguments to be passed to \code{\link{qr}}. }
}
\details{
  If \code{A} is square, the solution will be equivalent to (but should
  be faster than) \code{MPinv(A) \%*\% b}.

  No check is made on the validity of a specified \code{rank}
  argument.  If \code{rank} is left unspecified, the rank of
  \code{A} is determined numerically by \code{\link{qr}}.
  If the wrong rank is specified, failure (with an unpredictable
  error message) is almost certain.
}
\value{
  If \code{b} is a vector, the result is a vector of length n.  If
  \code{b} is a matrix, the result is a matrix with n rows and the
  same number of columns as \code{b}.

  In either case, the \code{"rank"} attribute of the result contains
  either the \code{rank} argument if specified, or the algebraic rank
  as calculated internally by \code{qr(A, ...)$rank}.
}
\references{ Lecture notes by James P Reilly, at
  http://www.ece.mcmaster.ca/~reilly/ee731/ch10.ps }
\author{ David Firth }

\seealso{ \code{\link{qr.solve}}, \code{\link{qr.coef}}, \code{\link{MPinv}} }
\examples{
set.seed(1)
x <- rnorm(3)
z <- rnorm(3)
xz <- x + z
X <- cbind(x, z, xz)
XX <- crossprod(X) # has rank 2
MPinv(XX, method = "chol")
qrSolve(XX)
b <- rnorm(3)
crossprod(MPinv(XX), b) 
qrSolve(XX, b)
crossprod(XX, qrSolve(XX, b))
fitted(lm(b ~ -1 + XX))
}
\keyword{ array }

