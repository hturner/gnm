\name{summary.gnm}
\alias{summary.gnm}
\alias{print.summary.gnm}
\title{ Summarize Generalized Nonlinear Model Fits }
\description{
\code{summary} method for objects of class \code{``gnm''}
}
\usage{
##  S3 method for class 'gnm'
\method{summary}{gnm}(object, dispersion = NULL, correlation = FALSE,
                      symbolic.cor = FALSE, ...)

## S3 method for class 'summary.gnm'
\method{print}{summary.gnm}(x, digits = max(3, getOption("digits") - 3),
                            symbolic.cor = x$symbolic.cor, ...)
}
\arguments{
  \item{object}{ an object of class \code{"gnm"}. }
  \item{x}{ an object of class \code{"summary.gnm"}. }
  \item{dispersion}{ the dispersion parameter for the fitting family. By
  default it is obtained from \code{object}. }
  \item{correlation}{ logical: if \code{TRUE}, the correlation matrix of
    the estimated parameters is returned. }
  \item{digits}{ the number of siginificant digits to use when
    printing. }
  \item{symbolic.cor}{ logical: if \code{TRUE}, the correlations are
    printed in a symbolic form rather than numbers (see \code{symnum}). }
  \item{\dots}{ further arguments passed to or from other methods. }
}
\details{
  \code{print.summary.gnm} prints the original call to \code{gnm}; a
  summary of the deviance residuals from the model fit; the coefficients
  of the model; the residual deviance; the Akaike's Information
  Criterion value; the number of main iterations performed, and if
  requested, the correlation matrix.
  
  Only the lower triangle of the correlation matrix is printed, to two
  decimal places; to see the full matrix print \code{summary(object,
  correlation = TRUE)$correlation} directly.
}
\value{
  \code{summary.gnm} returns an object of class \code{"summary.gnm"},
  which is a list with components
  \item{call }{ the \code{"call"} component from object.  }
  \item{terms }{ the \code{"terms"} component from object. }
  \item{family }{ the \code{"family"} component from object. }
  \item{deviance }{ the \code{"deviance"} component from object. }
  \item{aic }{ the \code{"aic"} component from object. }
  \item{df.residual }{ the \code{"df.residual"} component from object. }
  \item{iter }{ the \code{"iter"} component from object. }
  \item{deviance.resid }{ the deviance residuals, see \code{residuals.glm}. }
  \item{coefficients }{ the \code{"coefficients"} component from object.}
  \item{dispersion }{ either the supplied argument or the estimated
    dispersion if the latter is \code{NULL}. }
  \item{cov.unscaled }{ the unscaled (\code{dispersion = 1}) estimated
    covariance matrix of the estimated coefficients. }
  \item{cov.scaled }{ ditto, scaled by \code{dispersion}. }
  \item{correlation }{ (only if \code{correlation} is true) the
    estimated correlations of the estimated coefficients. }
  \item{symbolic.cor }{ (only if \code{correlation} is true) the value
    of the argument \code{symbolic.cor}. }
}
\author{ Heather Turner }
\seealso{ \code{\link{gnm}}, \code{\link{summary}}}
\examples{
##  Following on from example(gnm)
data(cautres)
set.seed(1)

##  Fit model as before
doubleUnidiff <- gnm(Freq ~ election:vote + election:class:religion +
                     Mult(Exp(election), religion:vote) +
                     Mult(Exp(election), class:vote), family = poisson,
                     data = cautres)

## Summarize results
summary(doubleUnidiff)
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
