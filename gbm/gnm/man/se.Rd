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
se(model, estimate = "all", checkEstimability = TRUE, ...)
}
\arguments{
  \item{model}{ a model object of class \code{"gnm"}}
  \item{estimate}{ specifies non-eliminated parameters or linear
    combinations of parameters for which to find standard errors. In the
    first case either \code{"all"}, a character vector of names, a
    numeric vector of indices or \code{"pick"} to select from a Tk
    dialog. In the second case coefficients given as a vector or the
    rows of a matrix, such that \code{NROW(estimate)} is equal to 
    \code{length(coef(model))} or \code{length(coef(model)) -
      model$eliminate}. }
  \item{checkEstimability}{ logical: should the estimability of
    all specified combinations be checked?}
  \item{\dots}{ possible further arguments for
    \code{\link{checkEstimable}} }
}
\note{
  In the case where \code{estimate} is a numeric vector, \code{se} will
  assume indices have been specified if all the values of
  \code{estimate} are in \code{1:(length(coef(model)) - model$eliminate)}.
}
\value{
A data frame with two columns:
  \item{estimate }{The estimated parameter combinations}
  \item{SE }{Their estimated standard errors}
and row names the same as the column names (if any) of \code{coefMatrix}.
}
\author{ David Firth }
\seealso{ \code{\link{gnm}}, \code{\link{getContrasts}},
  \code{\link{checkEstimable}}} 
\examples{
data(yaish)
set.seed(1)

## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(-1 + educ), -1 + orig:dest), family = poisson,
               data = yaish,  subset = (dest != 7))
## Deviance is 200.3

## Get estimate and se for the contrast between educ4 and educ5 in the
## UNIDIFF multiplier
educ4.pos <- grep("Mult.*educ4", names(coef(unidiff)))
mycontrast <- rep(0, length(coef(unidiff)))
mycontrast[educ4.pos] <- 1
mycontrast[educ4.pos + 1] <- -1
se(unidiff, mycontrast)

## Get all of the contrasts with educ5 in the UNIDIFF multipliers
getContrasts(unidiff, rev(grep("Mult.*educ", names(coef(unidiff)))))
}
\keyword{ models }
\keyword{ regression }
\keyword{nonlinear}
