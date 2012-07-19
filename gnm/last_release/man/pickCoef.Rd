\name{pickCoef}
\alias{pickCoef}
\title{ Get Indices or Values of Selected Model Coefficients  }
\description{
  Get the indices or values of a subset of non-eliminated coefficients selected
  via a Tk dialog or by pattern matching.
}
\usage{
pickCoef(object, pattern = NULL, value = FALSE, ...)
}
\arguments{
  \item{object}{ a model object. }
  \item{pattern}{ character string containing a regular expression or
    (with \code{fixed = TRUE}) a pattern to be matched exactly. If
    missing, a Tk dialog will open for coefficient selection.  }
  \item{value}{ if \code{FALSE}, a named vector of indices, otherwise
    the value of the selected coefficients. }
  \item{\dots}{ arguments to pass on to \link[relimp]{pickFrom} if
    \code{pattern} is missing, otherwise \code{grep}. In particular, \code{fixed
    = TRUE} specifies that \code{pattern} is a string to be matched as is.}
}
\value{
  If \code{value = FALSE} (the default), a named vector of indices,
  otherwise the values of the selected coefficients. If no coefficients
  are selected the returned value will be \code{NULL}.
}
\author{ Heather Turner }
\seealso{ \code{\link{regexp}}, \code{\link{grep}},
  \code{\link[relimp]{pickFrom}}, \code{\link{ofInterest}}}
\examples{
set.seed(1)

### Extract indices for use with ofInterest

data(yaish)
## fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(educ), orig:dest),
               family = poisson, data = yaish, subset = (dest != 7))

## set coefficients in first constituent multiplier as 'ofInterest'
## using regular expression
ofInterest(unidiff) <- pickCoef(unidiff, "[.]educ")

## summarise model, only showing coefficients of interest
summary(unidiff)

## get contrasts of these coefficients
getContrasts(unidiff, ofInterest(unidiff))

### Extract coefficients to use as starting values

## fit diagonal reference model with constant weights
set.seed(1)
data(voting)
## reconstruct counts voting Labour/non-Labour
count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)

classMobility <- gnm(yvar ~ -1 + Dref(origin, destination), 
                     family = binomial, data = voting)

## create factors indicating movement in and out of salariat (class 1)
upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

## extract diagonal effects from first model to use as starting values
diagCoef <- pickCoef(classMobility, "Dref(., .)", fixed = TRUE,
                     value = TRUE)

## fit separate weights for the "socially mobile" groups
## -- there are now 3 parameters for each weight
socialMobility <- gnm(yvar ~ -1 + Dref(origin, destination,
                                       delta = ~ 1 + downward + upward),
                      family = binomial, data = voting,
                      start = c(rep(NA, 6), diagCoef))
}
\keyword{ models }
