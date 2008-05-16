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
getContrasts(model, set = NULL,  refLevel = "first", dispersion = NULL,
  use.eliminate = TRUE, ...)
}
\arguments{
  \item{model}{ a model object of class \code{"gnm"}.}
  \item{set}{ a vector of indices (numeric) or coefficient names
    (character). If \code{NULL}, a Tk dialog will open for parameter
    selection. }
  \item{refLevel}{either a single numeric index, or a vector
    of real numbers which sum to 1, or one of the character
    strings \code{"first"}, \code{"last"} or \code{"mean"}.}
  \item{dispersion}{either \code{NULL}, or a positive number by which
    the model's variance-covariance matrix should be scaled.}
  \item{use.eliminate}{ logical; see \code{\link{vcov.gnm}}}
  \item{\dots}{ arguments to pass to other functions. }
}
\details{
  The indices in \code{set} must all be in \code{1:length(coef(object))}.  If
  \code{set = NULL}, a Tk dialog is presented for the selection
  of indices (model coefficients).

  For the set of coefficients selected, contrasts
  and their standard errors are computed.  
  A check is performed first on the estimability of all such
  contrasts.  The specific contrasts to be computed are controlled by the
  choice of \code{refLevel}: this may be \code{"first"} (the default),
  for contrasts with the first of the selected coefficients, or
  \code{"last"} for contrasts with the last, or \code{"mean"} for
  contrasts with the arithmetic mean of the coefficients in the selected
  set; or it may be an arbitrary vector of weights (summing
  to 1, not necessarily all non-negative)
  which specify a weighted mean against which contrasts are taken;
  or it may be a single index specifying one of the coefficients with
  which all contrasts should be taken.  Thus, for example,
  \code{refLevel = 1} is equivalent to \code{refLevel = "first"}, and
  \code{refLevel = c(1/3, 1/3, 1/3)} is equivalent to \code{refLevel =
    "mean"} when there are three coefficients in the selected \code{set}.

  The reported quasi-variances (and corresponding quasi standard errors), for
  details of which see Firth (2003) or Firth and Menezes (2004), are
  invariant to the choice of \code{refLevel}.
  
}
\value{
  An object of class
  \code{qv} --- see \code{\link[qvcalc]{qvcalc}}.
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
               Mult(Exp(educ), orig:dest),
               ofInterest = "[.]educ", family = poisson,
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
