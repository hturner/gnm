\name{Logit}
\alias{Logit}
\title{ Specify the Logit of a Predictor in a gnm Model Formula }
\description{
  A function of class \code{"nonlin"} to specify the logit
  of a predictor in the formula argument to \code{\link{gnm}}.
}
\usage{
Logit(expression, inst = NULL)
}
\arguments{
  \item{expression}{ a symbolic expression representing the (possibly
    nonlinear) predictor. }
  \item{inst}{ (optional) an integer specifying the instance number of
    the term. }
}
\details{
  The \code{expression} argument is interpreted as the right hand side
  of a formula in an object of class \code{"formula"}, except that an
  intercept term is not added by default. Any function of class
  \code{"nonlin"} may be used in addition to the usual operators and
  functions.  
}
\value{
  A list with the components required of a \code{"nonlin"} function:
  \item{ predictors }{the \code{expression} argument passed to \code{Logit}}
  \item{ term }{a function to create a deparsed mathematical expression
    of the term, given a label for the predictor.}
  \item{ call }{the call to use as a prefix for parameter labels. } 
}
\author{ Heather Turner }
\seealso{ \code{\link{gnm}}, \code{\link{formula}},
  \code{\link{nonlin.function}}}
\examples{
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
