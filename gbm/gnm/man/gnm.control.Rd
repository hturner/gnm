\name{gnm.control}
\alias{gnm.control}
\title{Set Control Parameters for \code{gnm}}
\description{
  Set parameters to control the \code{gnm} fitting algorithm.  
}
\usage{
gnm.control(epsilon = 1e-04, startit = 10, maxit = 500, trace = FALSE)
}
\arguments{
  \item{epsilon}{a positive numeric value specifying the tolerance level for
    convergence.}
  \item{startit}{a positive integer specifying the number of start-up iterations
    to perform.}
  \item{maxit}{ a positive integer specifying the maximum number of main
    iterations to perform.}
  \item{trace}{ a logical value indicating whether the deviance
    should be printed after each main iteration.}
}
\details{
  The \code{gnm} fitting algorithm consists of two stages.  In the start-up
  iterations, parameters are updated one-term at a time.  In the main
  iterations, all the parameters are updated in each iteration, until
  convergence is reached or the number or iterations reaches
  \code{maxit}.

  Convergence is judged by comparing the squared components of the score vector
  with corresponding elements of the diagonal of the Fisher information
  matrix: if the ratios are all less than \code{epsilon^2}, iterations cease.
}
\value{
  A list with the arguments as components.
}
\author{Heather Turner }
\seealso{ \code{\link{gnm}} }
\examples{
## A variation on example(gnm)...

data(cautres)

##  Fit a "double UNIDIFF" model with the religion-vote and class-vote
##  interactions both modulated by nonnegative election-specific
##  multipliers.

set.seed(1)
doubleUnidiff <- gnm(Freq ~ election:vote + election:class:religion
    + Mult(Exp(election), religion:vote) + Mult(Exp(election), class:vote),
    family = poisson, data = cautres, trace = TRUE, startit = 1)

## ...struggles with only one start-up iteration!
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
