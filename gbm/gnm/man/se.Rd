\name{se}
\alias{se}
\title{ Standard Errors of Linear Parameter Combinations in gnm Models }
\description{
  Computes approximate standard errors for (a selection of) individual
  parameters or one or more linear combinations of the parameters in a
  \code{\link{gnm}} (generalized nonlinear model) object.  By default, a
  check is made first on the estimability of each specified combination.
}
\usage{
se(model, estimate = ofInterest(model), checkEstimability = TRUE, Vcov =
NULL, dispersion = NULL, ...)
}
\arguments{
  \item{model}{ a model object of class \code{"gnm"}.}
  \item{estimate}{ (optional) specifies parameters or linear
    combinations of parameters for which to find standard errors. In the
    first case either a character vector of names, a
    numeric vector of indices or \code{"[?]"} to select from a Tk
    dialog. In the second case coefficients given as a vector or the
    rows of a matrix, such that \code{NROW(estimate)} is equal to 
    \code{length(coef(model))} or \code{length(coef(model)) -
      model$eliminate}. If missing, standard errors are returned for all
    parameters in the model.}
  \item{checkEstimability}{ logical: should the estimability of
    all specified combinations be checked?}
  \item{Vcov}{ either NULL, or a matrix }
  \item{dispersion}{ either NULL, or a positive number }
  \item{\dots}{ possible further arguments for
    \code{\link{checkEstimable}}. }
}
\note{
  In the case where \code{estimate} is a numeric vector, \code{se} will
  assume that indices have been specified if all the values of
  \code{estimate} are in \code{seq(length(coef(model))}.

  Where both \code{Vcov} and \code{dispersion} are supplied, the
  variance-covariance matrix of estimated model coefficients is taken to
  be \code{Vcov * dispersion}.
}
\value{
A data frame with two columns:
  \item{Estimate }{The estimated parameter combinations}
  \item{Std. Error }{Their estimated standard errors}
If available, the column names of \code{coefMatrix} will be used to name
the rows.
}
\author{ David Firth }
\seealso{ \code{\link{gnm}}, \code{\link{getContrasts}},
  \code{\link{checkEstimable}}, \code{\link{ofInterest}}} 
\examples{
data(yaish)
set.seed(1)

## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(educ), orig:dest),
               ofInterest = "[.]educ", family = poisson,
               data = yaish,  subset = (dest != 7))
## Deviance is 200.3

## Get estimate and se for the contrast between educ4 and educ5 in the
## UNIDIFF multiplier
mycontrast <- numeric(length(coef(unidiff)))
mycontrast[ofInterest(unidiff)[4:5]] <- c(1, -1)
se(unidiff, mycontrast)

## Get all of the contrasts with educ5 in the UNIDIFF multipliers
getContrasts(unidiff, rev(ofInterest(unidiff)))
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
