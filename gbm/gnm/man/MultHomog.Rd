\name{MultHomog}
\alias{MultHomog}
\title{``Plug-in'' Function to Fit Multiplicative Interactions
  with Homogeneous Effects}
\description{
  A function to provide the objects and functions required to fit
  multiplicative interactions with homogeneous effects in a generalized
  nonlinear model using \code{gnm}.  
}
\usage{
MultHomog(...)
}
\arguments{
  \item{\dots}{a comma-separated list of two or more expressions representing
    the factors in the interaction.}
}
\details{
  Designed to work as a plug-in function for \code{gnm}, \code{MultHomog}
  produces the objects required to fit a multiplicative interaction in which
  level effects are constrained to be equal across factors. If the given
  factors do not have exactly the same levels, the set of levels is taken to be
  the union of the factor levels, sorted into increasing order.
  
  To specify a homogeneous multiplicative interaction in the formula argument
  to \code{gnm}, the symbolic wrapper \code{Nonlin} must be used, with a call
  to \code{MultHomog} as the first argument. 
}
\value{
  A list with the components required of a \code{gnm} plug-in function:
  \item{labels}{a character vector of labels for the parameters in the
    interaction.}
  \item{predictor} {a function that takes estimates of the parameters in the
    interaction and returns the fitted values.}
  \item{localDesignFunction} {a function that takes a vector of estimates of
    the parameters in the interaction and returns
    the local design matrix.}
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
data(occupational.status)
gnm(Freq ~ origin + destination + Diag(origin, destination) +
            Nonlin(MultHomog(origin, destination)), family = poisson,
            data = occupational.status)
## Deviance is 32.56, 34 d.f.
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
