\name{voting}
\alias{voting}
\docType{data}
\title{Data on Social Mobility and the Labour Vote}
\description{
  Voting data from the 1987 British general election, cross-classified
  by the class of the head of household and the class of their father.
}
\usage{voting}
\format{
  A data frame with 25 observations on the following 4 variables.
  \describe{
    \item{\code{percentage}}{the percentage of the cell voting Labour.}
    \item{\code{total}}{the cell count.}
    \item{\code{origin}}{a factor describing the father's class with
      levels \code{1:5}.} 
    \item{\code{destination}}{a factor describing the head of
      household's class with levels \code{1:5}.} 
  }
}
\source{
  Clifford, P. and Heath, A. F. (1993) The Political Consequences of
  Social Mobility. \emph{J. Roy. Stat. Soc. A}, \bold{156(1)}, 51-61.
}
\examples{
### Examples from Clifford and Heath paper
### (Results differ slightly - possible transcription error in
### published data?)
set.seed(1)

## reconstruct counts voting Labour/non-Labour
count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)

## fit diagonal reference model with constant weights
classMobility <- gnm(yvar ~ -1 + Dref(origin, destination), 
                     family = binomial, data = voting)
DrefWeights(classMobility)

## create factors indicating movement in and out of salariat (class 1)
upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

## fit separate weights for the "socially mobile" groups
socialMobility <- gnm(yvar ~ -1 + Dref(origin, destination,
                                       delta = ~ 1 + downward + upward),
                      family = binomial, data = voting)
DrefWeights(socialMobility)

## fit separate weights for downwardly mobile groups only
downwardMobility <- gnm(yvar ~ -1 + Dref(origin, destination,
                                         delta = ~ 1 + downward),
                        family = binomial, data = voting)
DrefWeights(downwardMobility)
}
\keyword{datasets}
