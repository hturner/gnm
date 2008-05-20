\name{checkEstimable}
\alias{checkEstimable}
\title{ Check Whether One or More Parameter Combinations in a gnm Model
 are Identified }
\description{
  For each of a specified set of linear combinations of parameters from a
  \code{\link{gnm}} model, checks numerically whether the combination's
  estimate is invariant to re-parameterization of the model.
}
\usage{
checkEstimable(model, combMatrix = diag(seq(along = coef(model))),
                           tolerance = 1e6 * .Machine$double.eps)
}
\arguments{
  \item{model}{  a model object of class \code{"gnm"} }
  \item{combMatrix}{  numeric: either a vector of length the same as
    \code{length(coef(model))}, or a matrix with that number of rows.
    Coefficients of one or more linear combinations of the model's parameters.}
  \item{tolerance}{ numeric: a threshold value for detection of
    non-estimability.  }
}
\value{A logical vector of length equal to the number of parameter combinations
  tested; \code{NA} where a parameter combination is identically zero.}
\author{ David Firth }
\seealso{
  \code{\link{gnm}}, 
  \code{\link{se}},
  \code{\link{getContrasts}}
}
\references{
  Catchpole, E.A. and Morgan, B.J.T. (1997).  Detecting parameter
  redundancy. \emph{Biometrika}, \bold{84}, 187--196.
  }
\examples{
data(yaish)
set.seed(1)

## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(educ), orig:dest), family = poisson,
               data = yaish, subset = (dest != 7))

## Check whether multiplier contrast educ4 - educ5 is estimable
ofInterest(unidiff) <- pickCoef(unidiff, "[.]educ")
mycontrast <- numeric(length(coef(unidiff)))
mycontrast[ofInterest(unidiff)[4:5]] <- c(1, -1)
checkEstimable(unidiff, mycontrast)
## should be TRUE

## Check whether multiplier educ4 itself is estimable
mycontrast[ofInterest(unidiff)[5]] <- 0
checkEstimable(unidiff, mycontrast)
## should be FALSE -- only *differences* are identified here
}
\keyword{ models }
\keyword{ regression }
\keyword{nonlinear}
