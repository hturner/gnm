\name{Raise}
\alias{Raise}
\title{ Specify a Predictor Raised to a Constant Power in an gnm Model }
\description{
  A function of class \code{"nonlin"} to specify a predictor raised to a
  constant power in the formula argument to \code{\link{gnm}}.
}
\usage{
Raise(expression, power = 1, inst = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression}{  a symbolic expression representing the (possibly
    nonlinear) predictor. }
  \item{power}{ the power. }
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
  \item{ predictors }{the \code{expression} argument passed to \code{Raise}}
  \item{ term }{a function to create a deparsed mathematical expression
    of the term, given a label for the predictor.}
  \item{ call }{the call to use as a prefix for parameter labels. } 
}
\author{ Heather Turner }
\seealso{ \code{\link{gnm}}, \code{\link{formula}},
  \code{\link{nonlin.function}}, \code{\link{Inv}}}
\examples{
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
