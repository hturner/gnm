\name{Diag}
\alias{Diag}
\title{ Equality of two or more factors }
\description{
\code{Diag} is used mainly in regression models for data classified by 
two or more factors with the same levels.  By default, \code{Diag} 
operates on k-level factors to produce a new factor having k+1 levels;
if \code{binary = TRUE} is specified, the result is a coarser binary 
variable equal to 1 where all of the input factors are equal and 0 
otherwise.
}
\usage{
Diag(..., binary = FALSE)
}
\arguments{
  \item{\dots}{ One or more factors}
  \item{binary}{ Logical }
}
\value{
  Either a factor (if \code{binary = FALSE}) or a 0-1 numeric vector
  (if \code{binary = TRUE}).
}
\author{ David Firth }
\seealso{ \code{\link{Symm}}}
\examples{
row <- gl(4, 4, 16)
col <- gl(4, 1, 16)
Diag(row, col)
matrix(Diag(row, col, binary = TRUE), 4, 4)
}
\keyword{ models }
