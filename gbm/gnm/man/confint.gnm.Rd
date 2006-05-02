\name{confint.gnm}
\alias{confint.gnm}
\alias{confint.profile.gnm}
\title{ Compute Confidence Intervals of Parameters in a Generalized Nonlinear Model }
\description{
  Computes confidence intervals for one or more parameters in a
  generalized nonlinear model, based on the profiled deviance. 
}
\usage{
\method{confint}{gnm}(object, parm = ofInterest(object), level = 0.95,
    trace = FALSE, ...)

\method{confint}{profile.gnm}(object, parm = names(object), level = 0.95, ...)
}
\arguments{
  \item{object}{ an object of class \code{"gnm"} or \code{"profile.gnm"}}
  \item{parm}{ (optional) either a numeric vector of indices or a
    character vector of names, specifying the parameters for which
    confidence intervals are to be estimated. If \code{parm} is missing,
    confidence intervals are found for all parameters.} 
  \item{level}{ the confidence level required. }
  \item{trace}{ a logical value indicating whether profiling should be
    traced. }
  \item{\dots}{ arguments passed to or from other methods }
}
\details{
  These are methods for the generic function \code{confint} in the
  \code{base} package.

  For \code{"gnm"} objects, \code{profile.gnm} is first called to
  profile the deviance over each parameter specified by \code{parm}, or
  over all parameters in the model if \code{parm} is missing.

  The method for \code{"profile.gnm"} objects is then called, which
  interpolates the deviance profiles to estimate the limits of the
  confidence interval for each parameter, see \code{\link{profile.gnm}}
  for more details.

  If a \code{"profile.gnm"} object is passed directly to \code{confint},
  parametrs specified by \code{parm} must be a subset of the profiled
  parameters.
  
  For unidentified parameters a confidence interval cannot be calculated
  and the limits will be returned as \code{NA}. If the deviance curve
  has an asymptote and a limit of the confidence interval cannot be
  reached, the limit will be returned as \code{-Inf} or \code{Inf} 
  appropriately. If the range of the profile does not extend far enough
  to estimate a limit of the confidence interval, the limit will be
  returned as \code{NA}. In such cases, it may be desirable create a
  profile object directly, see \code{\link{profile.gnm}} for more
  details.
}
\value{
  A matrix (or vector) with columns giving lower and upper
  confidence limits for each parameter. These will be labelled as
  (1-level)/2 and 1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
}
\author{ Heather Turner }
\seealso{ \code{\link{profile.gnm}}, \code{\link{gnm}}, \code{\link{profile.glm}}}
\examples{
data(yaish)
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(-1 + educ), -1 + orig:dest), 
               ofInterest = "Mult.*educ", constrain = "Mult.*educ1",
               family = poisson, data = yaish, subset = (dest != 7))

## Letting 'confint' compute profile
confint(unidiff, trace = TRUE)
##                          2.5 \%     97.5 \%
## Mult1.Factor1.educ1         NA         NA
## Mult1.Factor1.educ2 -0.5978901  0.1022447
## Mult1.Factor1.educ3 -1.4836854 -0.2362378
## Mult1.Factor1.educ4 -2.5792398 -0.2953420
## Mult1.Factor1.educ5       -Inf -0.7007616

## Creating profile object first with user-specified stepsize
prof <- profile(unidiff, trace = TRUE, stepsize = 0.1)
confint(prof, ofInterest(unidiff)[2:5])
##                          2.5 \%     97.5 \%
## Mult1.Factor1.educ2 -0.5978324  0.1022441
## Mult1.Factor1.educ3 -1.4834753 -0.2362138
## Mult1.Factor1.educ4         NA -0.2950790
## Mult1.Factor1.educ5         NA         NA

## For 95\% confidence interval, need to estimate parameters for which
## z = +/- 1.96. Profile has not gone far enough for last two parameters
range(prof[[4]]$z)
range(prof[[5]]$z)
}
\keyword{ models }
\keyword{ nonlinear }