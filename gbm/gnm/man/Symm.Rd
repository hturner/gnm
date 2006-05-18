\name{Symm}
\alias{Symm}
\title{ Symmetric Interaction of Factors }
\description{
  \code{Symm} codes the symmetric interaction of factors having
  the same set of levels, for use in regression models of symmetry or 
  quasi-symmetry.
}
\usage{
Symm(...)
}
\arguments{
  \item{\dots}{ one or more factors. }
}
\value{
A factor whose levels index the symmetric interaction of all factors 
supplied as input.
}
\note{\code{Symm} relies on the gtools package from CRAN}
\author{ David Firth }
\seealso{ \code{\link{Diag}}}
\examples{
if (require(gtools)) {
    row <- gl(4, 4, 16)
    col <- gl(4, 1, 16)
    symm4by4 <- Symm(row, col)
    matrix(symm4by4, 4, 4)
} else cat("Symm requires the gtools package to be installed\n")
}
\keyword{ models }
