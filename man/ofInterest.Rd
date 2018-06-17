\name{ofInterest}
\alias{ofInterest}
\alias{ofInterest<-}
\title{ Coefficients of Interest in a Generalized Nonlinear Model }
\description{
  Retrieve or set the \code{"ofInterest"} component of a \code{"gnm"}
  (generalized nonlinear model) object.
}
\usage{
ofInterest(object)
ofInterest(object) <- value
}

\arguments{
  \item{object}{ an object of class \code{"gnm"}. }
  \item{value}{ a numeric vector of indices specifying the subset of
    (non-eliminated) coefficients of interest, or \code{NULL} to specify
    that all non-eliminated coefficients are of interest.  }  
}
\details{
  The \code{"ofInterest"} component of a \code{"gnm"} object is a named
  numeric vector of indices specifying a subset of the non-eliminated
  coefficients which are of specific interest.

  If the \code{"ofInterest"} component is non-NULL, printed summaries of
  the model only show the coefficients of interest. In addition 
  methods for \code{"gnm"} objects which may be applied to a subset of
  the parameters are by default applied to the coefficients of interest.

  These functions provide a way of extracting and replacing the
  \code{"ofInterest"} component. The replacement function prints the
  replacement value to show which parameters have been specified by
  \code{value}.  
}
\value{
  A named vector of indices, or \code{NULL}.
}
\author{ Heather Turner }
\note{
  Regular expression matching is performed using \code{grep} with
  default settings.
}
\seealso{ \code{\link{grep}}, \code{\link{gnm}}, \code{\link{se.gnm}},
  \code{\link{getContrasts}},\code{\link{profile.gnm}}, \code{\link{confint.gnm}}}
\examples{
set.seed(1)

## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(educ), orig:dest),
               ofInterest = "[.]educ", family = poisson,
               data = yaish,  subset = (dest != 7))
ofInterest(unidiff)

## Get all of the contrasts with educ1 in the UNIDIFF multipliers
getContrasts(unidiff, ofInterest(unidiff))

## Get estimate and se for the contrast between educ4 and educ5 in the
## UNIDIFF multiplier
mycontrast <- numeric(length(coef(unidiff)))
mycontrast[ofInterest(unidiff)[4:5]] <- c(1, -1)
se(unidiff, mycontrast)
}
\keyword{ models }
