\name{gnmConstrain}
\alias{gnmConstrain}
\title{ Re-run \code{gnm} Fit With Specified Constraints}
\description{
  Re-fit a generalized nonlinear model produced by \code{gnm},
  constraining specified parameters to zero.  
}
\usage{
gnmConstrain(model, constrain, start = coef(model), ...)
}
\arguments{
  \item{model}{ a \code{gnm} object. }
  \item{constrain}{ either a list of indices, or a logical vector,
    indicating elements of the parameter vector to set to zero. }
  \item{start}{ a vector of starting values for the parameters in the
    predictor. }
  \item{\dots}{ further arguments to pass to \code{gnm}. }
}
\value{
  an object of class \code{``gnm''}, see \code{gnm} for details.
}
\author{ Heather Turner }

\seealso{ \code{gnm} }
\examples{
##  Following on from example(gnm)...

data(cautres)

##  Fit model as before
doubleUnidiff <- gnm(Freq ~ election:vote + election:class:religion
    + Mult(Exp(election), religion:vote) + Mult(Exp(election), class:vote),
    family = poisson, data = cautres)

##  Re-fit model with Mult2.Factor1.election4 set to zero
doubleUnidiffConstrained <- gnmConstrain(doubleUnidiff, constrain = list(118))

##  Examine the multipliers of the class-vote log odds ratios

coef(doubleUnidiffConstrained)[118:121]

## ...as using 'getContrasts' in example(gnm) (to 4 d.p.).

}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
