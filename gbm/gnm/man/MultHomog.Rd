\name{MultHomog}
\alias{MultHomog}
\title{gnm Plug-in Function to Fit Multiplicative Interactions with
 Homogeneous Effects}
\description{
  A function to provide the objects and functions required to fit
  multiplicative interactions with homogeneous effects in a generalized
  nonlinear model using \code{\link{gnm}}.  
}
\usage{
MultHomog(...)
}
\arguments{
  \item{\dots}{ a comma-separated list of two or more factors. }
}
\details{
  Designed to work as a plug-in function for \code{gnm}, \code{MultHomog}
  produces the objects required to fit a multiplicative interaction with
  one component in which the constituent multipliers are the effects of
  two or more factors and the effects of these factors are constrained
  to be equal when the factor levels are equal. Thus the interaction
  effect would be
  \deqn{\gamma_i\gamma_j...}{gamma_i gamma_j ...}
  for an observation with level \eqn{i} of the first factor, level
  \eqn{j} of the second factor and so on, where \eqn{\gamma_l}{gamma_l}
  is the effect for level \eqn{l} of the homogeneous multiplicative factor.
  
  To specify a homogeneous multiplicative interaction in the formula argument
  to \code{gnm}, the symbolic wrapper \code{Nonlin} must be used, with a call
  to \code{MultHomog} as the first argument.

  If the factors passed to \code{MultHomog} do not have exactly the same
  levels, the set of levels is taken to be the union of the factor
  levels, sorted into increasing order. 
}
\value{
  A list with the components required of a \code{gnm} plug-in function:
  \item{labels}{ a character vector of labels for the parameters in the
    interaction. }
  \item{predictor}{ a function that takes estimates of the parameters in the
    interaction and returns the fitted values. }
  \item{localDesignFunction}{ a function that takes a vector of estimates of
    the parameters in the interaction and returns
    the local design matrix. }
}
\references{
  Goodman, L. A. (1979) Simple Models for the Analysis of Association in
  Cross-Classifications having Ordered
  Categories. \emph{J. Am. Stat. Assoc.}, \bold{74(367)}, 537-552. }
\author{ Heather Turner }
\seealso{\code{\link{gnm}} and \code{\link{Nonlin}} for usage.
  
  \code{\link{Dref}} for another \code{gnm} plug-in function.
  
  \code{\link{Mult}} for specifying multiplicative interactions in
  \code{gnm} formulae.}
\examples{
set.seed(1)
data(occupationalStatus)

##  Fit an association model with homogeneous row-column effects
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               Nonlin(MultHomog(origin, destination)), family = poisson,
               data = occupationalStatus)
## Deviance is 32.56, 34 d.f.
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
