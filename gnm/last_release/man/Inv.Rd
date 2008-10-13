\name{Inv}
\alias{Inv}
\title{ Specify the Reciprocal of a Predictor in a gnm Model Formula}
\description{
  A function of class \code{"nonlin"} to specify the reciprocal
  of a predictor in the formula argument to \code{\link{gnm}}.
}
\usage{
Inv(expression, inst = NULL)
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
  \item{ predictors }{the \code{expression} argument passed to \code{Inv}}
  \item{ term }{a function to create a deparsed mathematical expression
    of the term, given a label for the predictor.}
  \item{ call }{the call to use as a prefix for parameter labels. } 
}
\author{ Heather Turner }
\seealso{ \code{\link{gnm}}, \code{\link{formula}},
  \code{\link{nonlin.function}}}
\examples{
## One way to fit the logistic function without conditional
## linearity as in ?nls
library(gnm)
set.seed(1)
DNase1 <- subset(DNase, Run == 1)

test <- gnm(density ~ -1 +
            Mult(1, Inv(Const(1) + Exp(Mult(1 + offset(-log(conc)),
                                                Inv(1))))),
            start = c(NA, 0, 1), data = DNase1, trace = TRUE)
coef(test)
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }

