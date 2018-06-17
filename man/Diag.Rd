\name{Diag}
\alias{Diag}
\title{Equality of Two or More Factors}
\description{
Converts two or more factors into a new factor whose value is 0 where
the original factors are not all equal, and nonzero otherwise.
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
\details{
Used mainly in regression models for data classified by 
two or more factors with the same levels.  By default, 
operates on k-level factors to produce a new factor having k+1 levels;
if \code{binary = TRUE} is specified, the result is a coarser binary 
variable equal to 1 where all of the input factors are equal and 0 
otherwise.

If the original levels are identical the levels of the factor created 
in the \code{binary = FALSE} case will be in the same order, with 
\code{"."} added as the first level. Otherwise the levels of the new 
factor will be \code{"."} followed by the sorted combined levels.
}
\author{ David Firth and Heather Turner}
\seealso{\code{\link{Symm}}}
\examples{
rowfac <- gl(4, 4, 16)
colfac <- gl(4, 1, 16)
diag4by4 <- Diag(rowfac, colfac)
matrix(Diag(rowfac, colfac, binary = TRUE), 4, 4)
}
\keyword{ models }
