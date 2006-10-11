\name{voting}
\alias{voting}
\docType{data}
\title{Data on Social Mobility and the Labour Vote}
\description{
  Voting data from the 1987 British general election, cross-classified
  by the class of the head of household and the class of their father.
}
\usage{data(voting)}
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
data(voting)
## reconstruct counts voting Labour/non-Labour
count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)

## fit diagonal reference model with constant weights
classMobility <- gnm(yvar ~ -1 + Dref(origin, destination), 
                     family = binomial, data = voting)
prop.table(exp(coef(classMobility)[1:2]))

## create factors indicating movement in and out of salariat (class 1)
upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

## fit separate weights for the "socially mobile" groups
socialMobility <- gnm(yvar ~ -1 + Dref(origin, destination,
                                       delta = ~ 1 + downward + upward),
                      family = binomial, data = voting)
coefs <- coef(socialMobility)
baseline <- coefs[pickCoef(socialMobility, "(Intercept)")]
upwardTRUE <- coefs[pickCoef(socialMobility, "upwardTRUE")]
downwardTRUE <- coefs[pickCoef(socialMobility, "downwardTRUE")]
prop.table(exp(baseline))
prop.table(exp(upwardTRUE + baseline))
prop.table(exp(downwardTRUE + baseline))

## fit separate weights for downwardly mobile groups only
downwardMobility <- gnm(yvar ~ -1 + Dref(origin, destination,
                                         delta = ~ 1 + downward),
                        family = binomial, data = voting)
downwardMobility
coefs <- coef(downwardMobility)
baseline <- coefs[pickCoef(downwardMobility, "(Intercept)")]
downwardTRUE <- coefs[pickCoef(downwardMobility, "downwardTRUE")]
prop.table(exp(baseline))
prop.table(exp(downwardTRUE + baseline))
}
\keyword{datasets}
