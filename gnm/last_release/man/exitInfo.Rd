\name{exitInfo}
\alias{exitInfo}
\title{ Print Exit Information for gnm Fit }
\description{
  A utility function to print information on final iteration in
  \code{gnm} fit, intended for use when \code{gnm} has not converged.  
}
\usage{
exitInfo(object)
}
\arguments{
  \item{object}{ a \code{gnm} object. } 
}
\details{
  If \code{gnm} has not converged within the pre-specified maximum
  number of iterations, it may be because the algorithm has converged to
  a non-solution of the likelihood equations. In order to determine
  appropriate action, it is necessary to differentiate this case from
  one of near-convergence to the solution.

  \code{exitInfo} prints the absolute score and the corresponding
  convergence criterion for all parameters which failed to meet the
  convergence criterion at the last iteration. Clearly a small number of
  parameters with scores close to the criterion suggests
  near-convergence.  
}
\references{
  Vargas, M, Crossa, J, van Eeuwijk, F, Sayre, K D and Reynolds, M P
  (2001).  Interpreting treatment by environment interaction in agronomy
  trials.  \emph{Agronomy Journal} \bold{93}, 949--960.
}
\author{ Heather Turner }
\seealso{ \code{\link{gnm}}}
\examples{
data(cautres)

##  Fit a "double UNIDIFF" model with low iterMax for illustration!
set.seed(1)
doubleUnidiff <- gnm(Freq ~ election*vote + election*class*religion +
                     Mult(Exp(election), religion:vote) +
                     Mult(Exp(election), class:vote),
                     family = poisson, data = cautres, iterMax = 10)
exitInfo(doubleUnidiff)
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
