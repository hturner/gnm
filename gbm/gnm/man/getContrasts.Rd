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
  A list (of length \code{nSets}) of 2-component lists.  The first
  component, named \code{summary}, is a data frame containing variables
  \code{estimate}, \code{se} and \code{quasi.se}.  The \code{quasi.se}
  variable is
  not present if either the \code{qvcalc} package is unavailable
  or there are fewer than two estimable differences.  The second
  component, named \code{relative.errors}, is a character vector
  of length 2 containing the smallest and largest relative errors of the
  quasi-standard-error approximation, in the set of all simple
  contrasts; or \code{NULL} if quasi standard errors have not been
  calculated. See Firth (2003) or Firth and Menezes (2004)
  for details of quasi standard errors and their use.
  
}
\author{ David Firth }
\seealso{ \code{\link{gnm}}, \code{\link{se}},
  \code{\link{checkEstimable}}, \code{\link[qvcalc]{qvcalc}}}
\references{
  Firth, D (2003).  Overcoming the reference category problem in the
  presentation of statistical models.  \emph{Sociological Methodology}
  \bold{33}, 1--18.

  Firth, D and Menezes, R X de (2004).  Quasi-variances.
  \emph{Biometrika} \bold{91}, 65--80.
  }
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
