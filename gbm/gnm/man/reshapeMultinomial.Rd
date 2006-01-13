\name{reshapeMultinomial}
\alias{reshapeMultinomial}
\title{ Re-express Categorical Data as Counts }
\description{
  Expands the rows of a data frame containing a categorical variable
  \code{catvar}, such that \code{catvar} is replaced by a factor
  specifying the possible categories for each observation and a vector
  of 0/1 counts over these categories, whilst all other variables are
  replicated appropriately.
  %Expands the rows of a data frame containing a categorical variable
  %\code{catvar} with \eqn{c} possible categories, such that each
  %observation of \code{catvar} is represented by \eqn{c} 0/1 counts and
  %all other variables are replicated appropriately.
}
\usage{
reshapeMultinomial(data, catvar, countvar = "count", idvar = "id",
                   as.ordered = FALSE) 
}
\arguments{
  \item{data}{ a data frame. }
  \item{catvar}{ a character string specifying the factor in \code{data}
    which is to form the basis of the expansion. }
  \item{countvar}{ (optional) a character string to be used for the name
    of the new count variable }
  \item{idvar}{ (optional) a character string to be used for the name
    of the new factor identifying the original rows }
  \item{as.ordered}{ logical - whether the factor replacing
    \code{catvar} chould be of class \code{"ordered"} }
}
\details{
  Each row of the data frame is replicated \eqn{c} times, where \eqn{c}
  is the number of possible values that \code{catvar} can take. The
  expanded \code{catvar} factor is replaced by a factor of the same name
  specifying the \eqn{c} possible values for each case. The levels of
  the new factor will be in the same order as the original factor, but
  the new factor will not be of class \code{"ordered"} unless the
  argument \code{as.ordered = TRUE}. A new variable \code{countvar} is
  added which is equal to 1 for the observed category in each case and 0
  elsewhere. Finally a new factor \code{idvar} is added to index the
  cases.
}
\value{
  The expanded data frame as described in Details.
}
\author{ Heather Turner }
\note{
  Re-expressing categorical data in this way allows a multinomial
  response to be modelled as a poisson response, see examples.
}
\seealso{ \code{\link{gnm}}, \code{\link{multinom}},
  \code{\link{reshape}}, \code{\link[catspec]{mclgen}} }
\references{
  Anderson, J. A. (1984) Regression and Ordered Categorical
  Variables. \emph{J. R. Statist. Soc. B}, \bold{46(1)}, 1-30.
}
\examples{
### Example from ?multinom
data(birthwt, package = MASS)
summary(birthwt)
library(nnet)
bwt.mu <- multinom(low ~ . data = bwt)

## Equivalent using gnm - include unestimable main effects in model so 
## that interactions with low0 automatically set to zero, else could use
## 'constrain' argument. 
bwtLong <- reshapeMultinomial(bwt, "low")
bwt.po <- gnm(count ~  low*(.), eliminate = id, data = bwt2, family =
              "poisson") 
summary(bwt.po) # same deviance; df reflect extra id parameters

### Example from ?backPain
set.seed(1)
data(backPain)
summary(backPain)
backPainLong <- reshapeMultinomial(backPain, "pain")

## Fit models described in Table 5 of Anderson (1984)

noRelationship <- gnm(count ~ pain, eliminate = id,
                      family = "poisson", data = backPainLong)

oneDimensional <- update(noRelationship,
                         ~ . + Mult(pain - 1, x1 + x2 + x3 - 1))
}
\keyword{ manip }
\keyword{ models }
