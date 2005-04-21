\name{gnm.control}
\alias{gnm.control}
\title{Set Control Parameters for Model Fitting}
\description{
  Set parameters to control the \code{gnm} fitting algorithm.  
}
\usage{
gnm.control(epsilon = 1e-04, startit = 2, maxit = 500, trace = FALSE)
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
  iterations, any nonlinear parameters that are not specified by either the
  \code{start} argument of \code{gnm} or a plug-in function are
  updated one parameter at a time, then the linear parameters are
  jointly updated before the next iteration. In the main iterations, all
  the parameters are jointly updated, until convergence is reached or
  the number or iterations reaches \code{maxit}.

  Convergence is judged by comparing the squared components of the score vector
  with corresponding elements of the diagonal of the Fisher information
  matrix. If, for all components of the score vector, the ratio is less
  than \code{epsilon^2}, or the corresponding diagonal element of the
  Fisher information matrix is less than 1e-20, iterations cease.
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
system.time(
doubleUnidiff <- gnm(Freq ~ election:vote + election:class:religion
    + Mult(Exp(election - 1), religion:vote - 1) +
    Mult(Exp(election - 1), class:vote - 1),
    family = poisson, data = cautres, trace = TRUE))

## repeat with 10 start-up iterations
set.seed(1)
system.time(
doubleUnidiff <- gnm(Freq ~ election:vote + election:class:religion
    + Mult(Exp(election - 1), religion:vote - 1) +
    Mult(Exp(election - 1), class:vote - 1),
    family = poisson, data = cautres, trace = TRUE, startit = 10))
## ...at convergence by first main iteration, but much slower!
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
