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
### (Results differ slightly from paper - possible transcription error in
### published data?)
set.seed(1)
data(voting)
## reconstruct counts voting Labour/non-Labour
count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)

## fit diagonal reference model with constant weights
classMobility <- gnm(yvar ~ Nonlin(Dref(origin, destination)), 
                     family = binomial, data = voting)
prop.table(exp(coef(classMobility)[2:3]))
}

## create factors indicating movement in and out of salariat (class 1)
upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

## fit separate weights for the "socially mobile" groups
socialMobility <- gnm(yvar ~ Nonlin(Dref(origin, destination,
                                         formula = ~ 1 + downward + upward)),
                      family = binomial, data = voting)
prop.table(exp(coef(socialMobility)[c(4, 7)] + coef(socialMobility)[c(2, 5)]))
prop.table(exp(coef(socialMobility)[c(3, 6)] + coef(socialMobility)[c(2, 5)]))
prop.table(exp(coef(socialMobility)[c(2, 5)]))

## fit separate weights for downwardly mobile groups only
downwardMobility <- gnm(yvar ~ Nonlin(Dref(origin, destination,
                                           formula = ~ 1 + downward)),
                        family = binomial, data = voting)
downwardMobility
prop.table(exp(coef(downwardMobility)[c(3, 5)] + 
               coef(downwardMobility)[c(2, 4)]))
prop.table(exp(coef(downwardMobility)[c(2, 4)]))
}
\keyword{datasets}
