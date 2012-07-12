\name{gnm-defunct}
\alias{gnm-defunct}
\alias{Nonlin}
\alias{getModelFrame}
\alias{qrSolve}
\title{Defunct Functions in gnm Package}
\description{
  The functions listed here are no longer part of gnm as they are not
  needed any more. 
}
\usage{
Nonlin()
getModelFrame()
qrSolve(A, b, rank = NULL, ...)
}
\details{
  \code{Nonlin} is not needed any more as
  the plug-in architecture has been replaced by functions of class
  \code{"nonlin"}, see \code{\link{nonlin.function}}.

  \code{getModelFrame} was designed to work from within a plug-in
  function so is no longer needed.

  \code{qrSolve} was a function to solve the linear system Ax = b by two
  applications of QR decomposition. Alternative methods were found to be
  more robust.
}
\seealso{\code{\link{.Defunct}}}
\keyword{documentation}
