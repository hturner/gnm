\name{Symm}
\alias{Symm}
\title{ Symmetric Interaction of Factors }
\description{
  \code{Symm} codes the symmetric interaction of factors having
  the same set of levels, for use in regression models of symmetry or 
  quasi-symmetry.
}
\usage{
Symm(..., separator = ":")
}
\arguments{
  \item{\dots}{ one or more factors. }
  \item{separator}{ a character string of length 1 or more, to be used in naming the levels of the resulting interaction factor. }
}
\value{
A factor whose levels index the symmetric interaction of all factors 
supplied as input.
}
\author{ David Firth and Heather Turner }
\seealso{ \code{\link{Diag}}}
\examples{
# square table
rowfac <- gl(4, 4, 16)
colfac <- gl(4, 1, 16)
symm4by4 <- Symm(rowfac, colfac)
matrix(symm4by4, 4, 4)

# 3 x 3 x 3 table
ind <- expand.grid(A = 1:3, B = 1:3, C = 1:3)
symm3cubed <- with(ind, Symm(A, B, C))
array(symm3cubed, c(3, 3, 3))
}
\keyword{ models }
