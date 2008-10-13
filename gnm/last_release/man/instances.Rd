\name{instances}
\alias{instances}
\title{ Specify Multiple Instances of a Nonlinear Term in a gnm Model
  Formula } 
\description{
  A symbolic wrapper, for use in the formula argument to
  \code{\link{gnm}}, to specify multiple instances of a term specified
  by a function with an \code{inst} argument.
}
\usage{
instances(term, instances = 1)
}
\arguments{
  \item{term}{ a call to a function with an inst argument, which
    specifies some term. }
  \item{instances}{ the desired number of instances of the term. }
}
\value{
  A deparsed expression representing the summation of \code{term}
  specified with \code{inst = 1}, \code{inst = 2}, ..., \code{inst =
    instances}, which is used to create an expanded formula.
}
\author{ Heather Turner}
\seealso{\code{\link{gnm}}, \code{\link{formula}}, 
  \code{\link{nonlin.function}}, \code{\link{Mult}},
  \code{\link{MultHomog}}
}
\examples{
\dontrun{
## (this example can take quite a while to run)
##
## Fitting two instances of a multplicative interaction (i.e. a
## two-component interaction)
data(wheat)
yield.scaled <- wheat$yield * sqrt(3/1000)
treatment <- factor(paste(wheat$tillage, wheat$summerCrop, wheat$manure,
                          wheat$N, sep = ""))
bilinear2 <- gnm(yield.scaled ~ year + treatment +
                 instances(Mult(year, treatment), 2),
                 family = gaussian, data = wheat)
}
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
