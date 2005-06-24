\name{getContrasts}
\alias{getContrasts}
\title{ Estimated Contrasts and Standard Errors for Parameters in a gnm
 Model } 
\description{
For each set in a specified list of sets of parameters from a
\code{\link{gnm}} model, computes the estimated simple contrasts (i.e.,
differences) with the last parameter in the set, and estimated standard
errors for those estimated differences.
}
\usage{
getContrasts(model, sets = NULL, nSets = 1, ...)
}
\arguments{
  \item{model}{ a model object of class \code{"gnm"}}
  \item{sets}{ a vector of indices (if \code{nSets} is 1) or
    a list (of length \code{nSets}) of such vectors }
  \item{nSets}{ the number of vectors of indices to use}
  \item{\dots}{ arguments to pass to other functions }
}
\details{
  The indices must all be in \code{1:length(coef(object))}.  If
  \code{sets = NULL}, a Tk dialog is presented for the selection
  of indices (model coefficients).

  For each set of coefficients selected, differences with the last
  coefficient and their standard errors are computed.  
  A check is performed first on the estimability of all such
  differences.
}
\value{
  A list (of length \code{nSets}) of data frames, each containing variables
  \code{estimate} and \code{se}
}
\author{ David Firth }
\seealso{ \code{\link{gnm}}, \code{\link{se}}, \code{\link{checkEstimable}} }
\examples{
set.seed(1)
data(yaish)

## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ:orig + educ:dest +
               Mult(Exp(-1 + educ), orig:dest), family = poisson,
               data = yaish)
## Examine the education multipliers (differences on the log scale):
getContrasts(unidiff, grep("Mult1.Factor1", names(coef(unidiff))))
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
