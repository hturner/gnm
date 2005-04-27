\name{checkEstimable}
\alias{checkEstimable}
\title{ Check whether one or more parameter comninations in a gnm model
  is identified }
\description{
  For each of a specified set of parameter combinations, checks
  numerically whether the combination's estimate is invariant to
  re-parameterization of the model.
}
\usage{
checkEstimable(model, coefMatrix, tolerance = 1e-8)
}
\arguments{
  \item{model}{  a model object of class \code{"gnm"} }
  \item{coefMatrix}{  numeric: either a vector of length the same as
    \code{length(coef(model))}, or a matrix with that number of rows.
    Coefficients of one or more linear combinations of the model's parameters.}
  \item{tolerance}{ numeric: a threshold value for detection of
    non-estimability.  }
}
\value{A logical vector of length equal to the number of parameter combinations
  tested; \code{NA} where a parameter combination is identically zero.}
\references{
  Yaish, M (1998). Opportunities, Little Change. Class Mobility in
  Israeli Society: 1974-1991.  D.Phil. Thesis, Nuffield College,
  University of Oxford. 
}
\author{ David Firth }
\seealso{
  \code{\link{se}},
  \code{\link{getContrasts}}
}
\examples{
data(yaish)
set.seed(1)

## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ:orig + educ:dest +
               Mult(Exp(educ), orig:dest), family = poisson,
               data = yaish)

## Check whether Mult1.Factor1.educ4 - Mult1.Factor1.educ5 is estimable
educ4.pos <- grep("Mult1.Factor1.educ4", names(coef(unidiff)))
checkEstimable(unidiff, c(rep(0, educ4.pos - 1), 1, -1,
                          rep(0, length(coef(unidiff)) - educ4.pos - 1)))
## should be TRUE

## Check whether Mult1.Factor1.educ4 is estimable
checkEstimable(unidiff, c(rep(0, educ4.pos - 1), 1, 0,
                          rep(0, length(coef(unidiff)) - educ4.pos - 1)))
## should be FALSE -- only *differences* are identified here
}
\keyword{ models }
\keyword{ regression }
\keyword{nonlinear}
