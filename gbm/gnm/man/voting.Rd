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
    \item{\code{destination}}{a factor describing the son's class with
      levels \code{1:5}.} 
  }
}
\source{
  Clifford, P. and Heath, A. F. (1993) The Political Consequences of
  Social Mobility. \emph{J. Roy. Stat. Soc. A}, \bold{156(1)}, 51-61.
}
\examples{
data(voting)
}
\keyword{datasets}
