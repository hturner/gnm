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
checkEstimable(model, cmatrix, tolerance = 1e-8)
}
\arguments{
  \item{model}{  a model object of class \code{\dQuote{gnm}} }
  \item{cmatrix}{  numeric: either a vector of length the same as
    \code{length(coef(model))}, or a matrix with that number of rows.
    Coefficients of one or more linear combinations of the model's parameters.}
  \item{tolerance}{ numeric: a threshold value for detection of
    non-estimability.  }
}
\value{A logical vector of length equal to the number of parameter combinations
tested; \code{NA} where a parameter combination is identically zero.}
\author{ David Firth }
\seealso{
  \code{\link{se}},
  \code{\link{getContrasts}},
  \code{\link{gnmConstrain}}
}
\examples{
data(yaish)
set.seed(1)
yaish.unidiff <- gnm(Freq ~ educ:orig + educ:dest + Mult(Exp(educ),
    orig:dest), family = poisson, data = yaish)
educ4.pos <- grep("Mult1.Factor1.educ4", names(coef(yaish.unidiff)))
checkEstimable(yaish.unidiff,
 c(rep(0, educ4.pos - 1), 1, -1,
   rep(0, length(coef(yaish.unidiff)) - educ4.pos - 1)))
## should be TRUE
checkEstimable(yaish.unidiff,
 c(rep(0, educ4.pos - 1), 1, 0,
   rep(0, length(coef(yaish.unidiff)) - educ4.pos - 1)))
## should be FALSE -- only *differences* are identified here
}
\keyword{ models }
\keyword{ regression }
\keyword{nonlinear}
