\name{getContrasts}
\alias{getContrasts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{q}
\description{
q
}
\usage{
getContrasts(model, sets = NULL, nsets = 1,
             check.identifiability = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ a model object of class \code{"gnm"} }
  \item{sets}{ a vector of indices (if \code{nsets} is 1) or
    a list (of length nsets) of such vectors }
  \item{nsets}{ the number of vectors of indices to use}
  \item{check.identifiability}{logical: should
    identifiability of
    the chosen contrasts be checked first?}
  \item{\dots}{ arguments to pass to other functions }
}
\details{
  The indices must all be in \code{1:length(coef(object))}.  If
  \code{sets = NULL}, a Tk dialog is presented for the selection
  of indices (model coefficients).

  For each set of coefficients selected, differences with the last
  coefficient and their standard errors are computed.  By default
  check is preformed first on the identifiability of all such
  differences.
}
\value{
  A list (of length \code{nsets}) of data frames, each containing variables
  \code{estimate} and \code{se}
}
\author{ David Firth }
\seealso{ \code{\link{se}}, \code{\link{gnmConstrain}} }
\examples{
data(yaish)
set.seed(1)
yaish.unidiff <- gnm(Freq ~ educ:orig + educ:dest + Mult(Exp(educ),
    orig:dest), family = poisson, data = yaish)
## Examine the education multipliers (differences on the log scale):
getContrasts(yaish.unidiff, grep("Mult1.Factor1", names(coef(yaish.unidiff))))
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
