\name{se}
\alias{se}
\title{ Standard errors of linear parameter combinations in gnm models }
\description{
  Computes approximate standard error for one or more linear
  combinations of the parameters in a \code{\link{gnm}} (generalized
    nonlinear model) object.  By default, a check is made first on the
    estimability of each specified combination.
}
\usage{
se(model, cmatrix, check.estimability = TRUE, ...)
}
\arguments{
  \item{model}{  a model object of class \code{\dQuote{gnm}}}
  \item{cmatrix}{ numeric: either a vector of length the same as
    \code{length(coef(model))}, or a matrix with that number of rows.
    Coefficients of one or more linear combinations of the model's parameters.}
  \item{check.estimability}{ logical: should the estimability of
    all specified combinations be checked?}
  \item{\dots}{ possible further arguments for
    \code{\link{checkEstimable}} }
}
\value{
A data frame with two columns:
  \item{estimate }{The estimated parameter combinations}
  \item{se }{Their estimated standard errors}
and row names the same as the column names (if any) of \code{cmatrix}.
}
\author{ David Firth }
\seealso{ \code{\link{getContrasts}}, \code{\link{gnmConstrain}} }
\examples{
data(yaish)
set.seed(1)
yaish.unidiff <- gnm(Freq ~ educ:orig + educ:dest + Mult(Exp(-1 + educ),
    orig:dest), family = poisson, data = yaish)
## Deviance is 208.3
##
## Get estimate and se for the contrast between educ4 and educ5 in the
## UNIDIFF multiplier
##
educ4.pos <- grep("Mult.*educ4", names(coef(yaish.unidiff)))
mycontrast <- rep(0, length(coef(yaish.unidiff)))
mycontrast[educ4.pos] <- 1
mycontrast[educ4.pos + 1] <- -1
se(yaish.unidiff, mycontrast)
##
## Get all of the contrasts with educ5 in the UNIDIFF multipliers
##
getContrasts(yaish.unidiff, grep("Mult.*educ", names(coef(yaish.unidiff))))
}
\keyword{ models }
\keyword{ regression }
\keyword{nonlinear}
