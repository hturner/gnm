\name{residSVD}
\alias{residSVD}
\title{ Multiplicative Approximation of Model Residuals }
\description{
  This function uses the first \code{d} components of the singular value
  decomposition in order to approximate a vector of model residuals by a
  sum of \code{d} multiplicative terms, with the multiplicative
  structure determined by two specified factors.  It applies to models
  of class \code{lm}, \code{glm} or \code{gnm}.
}
\usage{
residSVD(model, fac1, fac2, d = 1)
}
\arguments{
  \item{model}{ an object of class \code{gnm}, \code{glm} or \code{lm} }
  \item{fac1}{ a factor }
  \item{fac2}{ a factor }
  \item{d}{ integer, the number of multiplicative terms to use in the approximation }
}
\details{
  This function operates on the matrix of mean residuals, with rows
  indexed by \code{fac1} and columns indexed by \code{fac2}.  For
  \code{glm} and \code{glm} models, the matrix entries are weighted
  working residuals.  The primary use of \code{residSVD} is to
  generate good starting values for the parameters in \code{\link{Mult}} terms
  in models to be fitted using \code{\link{gnm}}.   
}
\value{
  If \code{d = 1}, a numeric vector; otherwise a numeric
  matrix with \code{d} columns.
}
\author{ David Firth and Heather Turner }
\seealso{ \code{\link{gnm}}, \code{\link{Mult}}} 
\examples{
set.seed(1)
data(mentalHealth)
##  Goodman RC1  association model fits well (deviance 3.57, df 8)
mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)
## independence model
indep <- gnm(count ~ SES + MHS, family = poisson, data = mentalHealth)
mult1 <- residSVD(indep, SES, MHS)
## Now use mult1 as starting values for the RC1 association parameters
RC1model <- update(indep, . ~ . + Mult(SES, MHS),
                   start = c(coef(indep), mult1), trace = TRUE)
##  Similarly for the RC2 model:
mult2 <- residSVD(indep, SES, MHS, d = 2)
RC2model <- update(indep, . ~ . + instances(Mult(SES, MHS), 2),
                   start = c(coef(indep), mult2), trace = TRUE)
##
## See also example(House2001), where good starting values matter much more!
##
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
