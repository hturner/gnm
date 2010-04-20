\name{gnm-deprecated}
\alias{gnm-deprecated}
\alias{qrSolve}
\title{Deprecated Functions in gnm Package}
\description{
  These functions are provided for compatibility with older versions
  of the gnm package only, and may be defunct as soon as the next release.
}
\usage{
getModelFrame()
qrSolve(A, b, rank = NULL, ...)
}
\arguments{
  \item{A}{ A numeric matrix (m by n, say) }
  \item{b}{ A numeric vector of length m, or a matrix with m rows }
  \item{rank}{ Numeric: the algebraic rank of \code{A}, if known. }
  \item{\dots}{ Other possible arguments to be passed to \code{\link{qr}}. }
}
\details{
  The original help pages are available at
  \code{help("oldName-deprecated")} (note the quotes). 
}
\seealso{\code{\link{.Deprecated}}}
\keyword{documentation}
