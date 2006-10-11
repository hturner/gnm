\name{MultHomog}
\alias{MultHomog}
\title{Specify a Multiplicative Interaction with Homogeneous Effects in
  a gnm Model Formula}
\description{
  A function of class \code{"nonlin"} to specify a multiplicative
  interaction with homogeneous effects in the formula argument to
  \code{\link{gnm}}.   
}
\usage{
MultHomog(..., inst = NULL)
}
\arguments{
  \item{\dots}{ a comma-separated list of two or more factors. }
  \item{inst}{ (optional) an integer specifying the instance number of
    the term. }
}
\details{
  \code{MultHomog} specifies instances of a multiplicative interaction
  in which the constituent multipliers are the effects of two or more
  factors and the effects of these factors are constrained to be equal
  when the factor levels are equal.  Thus the interaction effect would
  be
  \deqn{\gamma_i\gamma_j...}{gamma_i gamma_j ...}
  for an observation with level \eqn{i} of the first factor, level
  \eqn{j} of the second factor and so on, where \eqn{\gamma_l}{gamma_l}
  is the effect for level \eqn{l} of the homogeneous multiplicative
  factor. 

  If the factors passed to \code{MultHomog} do not have exactly the same
  levels, the set of levels is taken to be the union of the factor
  levels, sorted into increasing order. 
}
\value{
   A list with the anticipated components of a \code{"nonlin"} function:
   \item{ predictors }{ the factors passed to \code{MultHomog}}
   \item{ common }{ an index to specify that common effects are to be
     estimated across the factors }
  \item{ term }{ a function to create a deparsed mathematical expression
    of the term, given labels for the predictors.}
  \item{ call }{ the call to use as a prefix for parameter labels. } 
}
\references{
  Goodman, L. A. (1979) Simple Models for the Analysis of Association in
  Cross-Classifications having Ordered
  Categories. \emph{J. Am. Stat. Assoc.}, \bold{74(367)}, 537-552. }
\note{Currently, \code{MultHomog} can only be used to specify a
  one-dimensional interaction. See examples for a workaround to specify
  interactions with more than one dimension. }
\author{ Heather Turner }
\seealso{\code{\link{gnm}}, \code{\link{formula}}, \code{\link{instances}},
  \code{\link{nonlin.function}}, \code{\link{Mult}}}
\examples{
set.seed(1)
data(friend)

###  Fit an association model with homogeneous row-column effects
rc1 <- gnm(Freq ~ r + c + Diag(r,c) + MultHomog(r, c),
           family = poisson, data = friend)
rc1

\dontrun{
###  (Takes a little while to run)
###  Extend to two-component interaction
rc2 <- update(rc1, . ~ . + MultHomog(r, c, inst = 2))
rc2
}
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
