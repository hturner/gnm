\name{expandCategorical}
\alias{expandCategorical}
\title{ Expand Data Frame by Re-expressing Categorical Data as Counts }
\description{
  Expands the rows of a data frame by re-expressing observations of a
  categorical variable specified by \code{catvar}, such that the
  column(s) corresponding to \code{catvar} are replaced by a factor 
  specifying the possible categories for each observation and a vector
  of 0/1 counts over these categories.
  %Expands the rows of a data frame containing a categorical variable
  %\code{catvar} with \eqn{c} possible categories, such that each
  %observation of \code{catvar} is represented by \eqn{c} 0/1 counts and
  %all other variables are replicated appropriately.
}
\usage{
expandCategorical(data, catvar, sep = ".", countvar = "count",
                  idvar = "id", as.ordered = FALSE, group = TRUE) 
}
\arguments{
  \item{data}{ a data frame. }
  \item{catvar}{ a character vector specifying factors in \code{data}
    whose interaction will form the basis of the expansion. }
  \item{sep}{ a character string used to separate the concatenated
    values of \code{catvar} in the name of the new interaction factor. } 
  \item{countvar}{ (optional) a character string to be used for the name
    of the new count variable. }
  \item{idvar}{ (optional) a character string to be used for the name
    of the new factor identifying the original rows (cases). }
  \item{as.ordered}{ logical - whether the new interaction factor should
    be of class \code{"ordered"}.}
  \item{group}{logical: whether or not to group individuals with common values
    over all covariates. } 
}
\details{
  Each row of the data frame is replicated \eqn{c} times, where \eqn{c}
  is the number of levels of the interaction of the factors specified by
  \code{catvar}. In the expanded data frame, the columns specified by
  \code{catvar} are replaced by a factor specifying the \eqn{r} possible
  categories for each case, named by the concatenated values of
  \code{catvar} separated by \code{sep}. The ordering of factor levels
  will be preserved in the creation of the new factor, but this factor
  will not be of class \code{"ordered"} unless the argument
  \code{as.ordered = TRUE}. A variable with name \code{countvar} is added
  to the data frame which is equal to 1 for the observed category in each
  case and 0 elsewhere. Finally a factor with name \code{idvar} is added
  to index the cases. }
\value{
  The expanded data frame as described in Details.
}
\author{ Heather Turner }
\note{
  Re-expressing categorical data in this way allows a multinomial
  response to be modelled as a poisson response, see examples.
}
\seealso{ \code{\link{gnm}}, \code{\link[nnet]{multinom}},
  \code{\link{reshape}}, \code{\link[catspec:mcl]{mclgen}} }
\references{
  Anderson, J. A. (1984) Regression and Ordered Categorical
  Variables. \emph{J. R. Statist. Soc. B}, \bold{46(1)}, 1-30.
}
\examples{
### Example from help(multinom, package = "nnet")
library(MASS)
example(birthwt)
library(nnet)
bwt.mu <- multinom(low ~ ., data = bwt)

## Equivalent using gnm - include unestimable main effects in model so 
## that interactions with low0 automatically set to zero, else could use
## 'constrain' argument. 
bwtLong <- expandCategorical(bwt, "low", group = FALSE)
bwt.po <- gnm(count ~  low*(. - id), eliminate = id, data = bwtLong, family =
              "poisson") 
summary(bwt.po) # same deviance; df reflect extra id parameters

### Example from ?backPain
set.seed(1)
summary(backPain)
backPainLong <- expandCategorical(backPain, "pain")

## Fit models described in Table 5 of Anderson (1984)

noRelationship <- gnm(count ~ pain, eliminate = id,
                      family = "poisson", data = backPainLong)

oneDimensional <- update(noRelationship,
                         ~ . + Mult(pain, x1 + x2 + x3))
}
\keyword{ manip }
\keyword{ models }
