\name{pickCoef}
\alias{pickCoef}
\title{ Get Indices For Selected Model Coefficients  }
\description{
  Get the indices for a subset of model coefficients selected via a Tk
  dialog or by regular expression matching.
}
\usage{
pickCoef(object, regexpr = NULL, ...)
}
\arguments{
  \item{object}{ a model object. }
  \item{regexpr}{ character string containing a regular expression. If
    missing, a Tk dialog will open for coefficient selection.  }
  \item{\dots}{ arguments to pass on to relimp:::pickFrom. }
}
\value{
  A named vector of indices, or \code{NULL}.
}
\author{ Heather Turner }
\note{
  Regular expression matching is performed using \code{grep} with
  default settings.
}
\seealso{ \code{\link{regexp}}, \code{\link{grep}}, \code{\link{ofInterest}}}
\examples{
set.seed(1)

### Use with ofInterest
data(yaish)
## fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(educ), orig:dest),
               family = poisson, data = yaish, subset = (dest != 7))

## set coefficients in first constituent multiplier as 'ofInterest'
ofInterest(unidiff) <- pickCoef(unidiff, "Mult.*educ")

## summarise model, only showing coefficients of interest
summary(unidiff)

## get contrasts of these coefficients
getContrasts(unidiff, ofInterest(unidiff))

### Use to pick out coefficients relating to different groups
data(voting)
## reconstruct counts voting Labour/non-Labour
count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)
## create factors indicating movement in and out of salariat (class 1)
upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

## fit separate weights for the "socially mobile" groups
socialMobility <- gnm(yvar ~ -1 + Dref(origin, destination,
                                       delta = ~ 1 + downward + upward),
                      family = binomial, data = voting)
## calculate weights
DrefWeights(socialMobility)
}
\keyword{ models }
