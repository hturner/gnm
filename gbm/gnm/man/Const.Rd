\name{Const}
\alias{Const}
\title{ Specify a Constant in a gnm Model Formula }
\description{
  A symbolic wrapper to specify a constant in the formula argument to
  \code{gnm}
}
\usage{
Const(const)
}
\arguments{
  \item{const}{ a numeric value. }
}
\details{
  \code{Const} is equivalent to \code{offset(rep(const, nObs))} where
  \code{nObs} is the number of observations being modelled.
}
\value{
  A call to \code{rep} used to create a variable representing the
  constant in the model frame.
}
\author{ Heather Turner }
\seealso{\code{\link{gnm}}, \code{\link{formula}}, \code{\link{offset}}}
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
