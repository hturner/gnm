\name{predict.gnm}
\alias{predict.gnm}
\title{ Predict Method for Generalized Nonlinear Models }
\description{
  Obtains predictions and optionally estimates standard errors of
  those predictions from a fitted generalized nonlinear model object.
}
\usage{
\method{predict}{gnm}(object, newdata = NULL,
type = c("link", "response", "terms"), se.fit = FALSE, dispersion =
NULL, terms = NULL, na.action = na.exclude, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ a fitted object of class inheriting from \code{"gnm"}. }
  \item{newdata}{ optionally, a data frame in which to look for variables with
          which to predict.  If omitted, the fitted predictors are used. }
  \item{type}{ the type of prediction required.  The default is on the scale
          of the predictors; the alternative \code{"response"} is on
          the scale of the response variable.  Thus for a default
          binomial model the default predictions are of log-odds
          (probabilities on logit scale) and \code{type = "response"} gives
          the predicted probabilities.  The \code{"terms"} option returns a
          matrix giving the fitted values of each term in the model
          formula on the predictor scale.

	  The value of this argument can be abbreviated. }
  \item{se.fit}{ logical switch indicating if standard errors are required. }
  \item{dispersion}{  the dispersion of the fit to be assumed in computing
          the standard errors.  If omitted, that returned by \code{summary}
          applied to the object is used. }
  \item{terms}{ with \code{type="terms"} by default all terms are returned. A
          character vector specifies which terms are to be returned }
  \item{na.action}{ function determining what should be done with missing values
          in \code{newdata}.  The default is to predict \code{NA}. }
  \item{\dots}{ further arguments passed to or from other methods. }
}
\details{
  If \code{newdata} is omitted the predictions are based on the data used
  for the fit.  In that case how cases with missing values in the
  original fit is determined by the \code{na.action} argument of that
  fit.  If \code{na.action = na.omit} omitted cases will not appear in
  the residuals, whereas if \code{na.action = na.exclude} they will
  appear (in predictions and standard errors), with residual value
  \code{NA}.  See also \code{\link{napredict}}.
}
\value{
  If \code{se = FALSE}, a vector or matrix of predictions.  If \code{se =
  TRUE}, a list with components
  \item{ fit }{ predictions.}
  \item{ se.fit }{ estimated standard errors.}
  \item{ residual.scale }{ a scalar giving the square root of the dispersion used
          in computing the standard errors.}
}
\references{ Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S }}
\author{ Heather Turner }
\note{Variables are first looked for in 'newdata' and then searched for
     in the usual way (which will include the environment of the
     formula used in the fit).  A warning will be given if the
     variables found are not of the same length as those in 'newdata'
     if it was supplied.}
\seealso{ \code{\link{gnm}} }
\examples{
set.seed(1)

##  Fit an association model with homogeneous row-column effects
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), family = poisson,
               data = occupationalStatus)

## Fitted values (expected counts)
predict(RChomog, type = "response", se.fit = TRUE)

## Fitted values on log scale
predict(RChomog, type = "link", se.fit = TRUE)
}

\keyword{ models }
\keyword{ nonlinear }
