\name{Symm}
\alias{Symm}
\title{ Symmetric interaction of factors }
\description{
  \code{Symm} codes the symmetric interaction of factors having
  the same set of levels, for use in regression models of symmetry or 
  quasi-symmetry.
}
\usage{
Symm(...)
}
\arguments{
  \item{\dots}{ One or more factors }
}
\value{
A factor whose levels index the symmetric interaction of all factors 
supplied as input.
}
\author{ David Firth }
\seealso{ \code{\link{Diag}}}
\examples{
row <- gl(4, 4, 16)
col <- gl(4, 1, 16)
\dontrun{
    Symm(row, col)
    matrix(Symm(row, col), 4, 4)
}
}
\keyword{ models }