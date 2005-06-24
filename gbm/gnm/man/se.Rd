\name{se}
\alias{se}
\title{ Standard Errors of Linear Parameter Combinations in gnm Models }
\description{
  Computes approximate standard error for one or more linear
  combinations of the parameters in a \code{\link{gnm}} (generalized
  nonlinear model) object.  By default, a check is made first on the
  estimability of each specified combination.
}
\usage{
se(model, coefMatrix, checkEstimability = TRUE, ...)
}
\arguments{
  \item{model}{  a model object of class \code{"gnm"}}
  \item{coefMatrix}{ numeric: either a vector of length the same as
    \code{length(coef(model))}, or a matrix with that number of rows.
    Coefficients of one or more linear combinations of the model's parameters.}
  \item{checkEstimability}{ logical: should the estimability of
    all specified combinations be checked?}
  \item{\dots}{ possible further arguments for
    \code{\link{checkEstimable}} }
}
\value{
A data frame with two columns:
  \item{estimate }{The estimated parameter combinations}
  \item{se }{Their estimated standard errors}
and row names the same as the column names (if any) of \code{coefMatrix}.
}
\author{ David Firth }
\seealso{ \code{\link{gnm}}, \code{\link{getContrasts}},
  \code{\link{checkEstimable}}} 
\examples{
data(yaish)
set.seed(1)

## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ:orig + educ:dest +
               Mult(Exp(-1 + educ), orig:dest), family = poisson,
               data = yaish)
## Deviance is 208.3

## Get estimate and se for the contrast between educ4 and educ5 in the
## UNIDIFF multiplier
educ4.pos <- grep("Mult.*educ4", names(coef(unidiff)))
mycontrast <- rep(0, length(coef(unidiff)))
mycontrast[educ4.pos] <- 1
mycontrast[educ4.pos + 1] <- -1
se(unidiff, mycontrast)

## Get all of the contrasts with educ5 in the UNIDIFF multipliers
getContrasts(unidiff, grep("Mult.*educ", names(coef(unidiff))))
}
\keyword{ models }
\keyword{ regression }
\keyword{nonlinear}
