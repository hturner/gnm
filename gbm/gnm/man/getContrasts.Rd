\name{getContrasts}
\alias{getContrasts}
\title{ Estimated Contrasts and Standard Errors for Parameters in a gnm
 Model } 
\description{
For each set in a specified list of sets of parameters from a
\code{\link{gnm}} model, computes the estimated simple contrasts (i.e.,
differences) with the first parameter in the set, and estimated standard
errors for those estimated differences.  Where possible, quasi standard
errors are also computed.
}
\usage{
getContrasts(model, sets = NULL, nSets = 1, dispersion = NULL,
  use.eliminate = TRUE, ...)
}
\arguments{
  \item{model}{ a model object of class \code{"gnm"}.}
  \item{sets}{ a vector of indices or
    a list of such vectors. If \code{NULL}, a Tk
    dialog will open for parameter selection. }
  \item{nSets}{ the number of vectors of indices to use when
    \code{sets} is \code{NULL}.}
  \item{dispersion}{either \code{NULL}, or a positive number by which
    the model's variance-covariance matrix should be scaled.}
  \item{use.eliminate}{ logical; see \code{\link{vcov.gnm}}}
  \item{\dots}{ arguments to pass to other functions. }
}
\details{
  The indices must all be in \code{1:length(coef(object))}.  If
  \code{sets = NULL}, a Tk dialog is presented for the selection
  of indices (model coefficients).

  For each set of coefficients selected, differences with the first
  coefficient and their standard errors are computed.  
  A check is performed first on the estimability of all such
  differences.

  If \code{sets} is non-\code{NULL}, the value of the \code{nsets}
  argument is ignored.
}
\value{
  EITHER a list (normally of length \code{nSets}) of objects of class
  \code{qv} --- see \code{\link[qvcalc]{qvcalc}} --- (when \code{sets}
  is a list, or when \code{sets} is \code{NULL} and \code{nSets > 1});

  OR an object of class \code{qv} (otherwise).
}
\author{ David Firth }
\seealso{ \code{\link{gnm}}, \code{\link{se}},
  \code{\link{checkEstimable}}, \code{\link[qvcalc]{qvcalc}},
  \code{\link{ofInterest}}}
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

## Fit the "UNIDIFF" mobility model across education levels -- see ?yaish
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(-1 + educ), -1 + orig:dest),
               ofInterest = "Mult.*educ", family = poisson,
               data = yaish,  subset = (dest != 7))
## Examine the education multipliers (differences on the log scale):
unidiffContrasts <- getContrasts(unidiff, ofInterest(unidiff))
plot(unidiffContrasts,
  main = "Unidiff multipliers (log scale): intervals based on
           quasi standard errors",
  xlab = "Education level", levelNames = 1:5)
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
