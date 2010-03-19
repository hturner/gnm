\name{parameters}
\alias{parameters}
\title{ Extract Constrained and Estimated Parameters from a gnm Object}
\description{
  A function to extract non-eliminated parameters from a \code{"gnm"}
  object, including parameters that were constrained.
}
\usage{
parameters(object)
}
\arguments{
  \item{object}{ an object of class \code{"gnm"}. }
}
\details{
  \code{parameters} acts like \code{coefficients} except that for
  constrained parameters, the value at which the parameter was
  constrained is returned instead of \code{NA}.
}
\value{
  A vector of parameters.
}
\author{ Heather Turner }
\seealso{ \code{\link{coefficients}}, \code{\link{gnm}} }
\examples{
data(occupationalStatus)
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), family = poisson,
               data = occupationalStatus, ofInterest = "MultHomog",
               constrain = "MultHomog.*1")
coefficients(RChomog)
parameters(RChomog)
}
\keyword{ models }
