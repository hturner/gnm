\name{summary.gnm}
\alias{summary.gnm}
\alias{print.summary.gnm}
\title{ Summarize Generalized Nonlinear Model Fits }
\description{
\code{summary} method for objects of class \code{"gnm"}
}
\usage{
\method{summary}{gnm}(object, dispersion = NULL, correlation = FALSE,
                      symbolic.cor = FALSE, with.eliminate = FALSE, ...)

\method{print}{summary.gnm}(x, digits = max(3, getOption("digits") - 3),
                            signif.stars = getOption("show.signif.stars"),
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
    printed in a symbolic form rather than numbers (see
    \code{symnum}). }
  \item{signif.stars}{ logical. If \code{TRUE}, "significance stars" are
    printed for each coefficient. }
  \item{with.eliminate}{ Logical. If \code{TRUE}, any eliminated
    coefficients are included in the summary. }
  \item{\dots}{ further arguments passed to or from other methods. }
}
\details{
  \code{print.summary.gnm} prints the original call to \code{gnm}; a
  summary of the deviance residuals from the model fit; the coefficients
  of the model; the residual deviance; the Akaike's Information
  Criterion value, and the number of main iterations performed.

  Standard errors, z-values and p-values are printed alongside the
  coefficients, with "significance stars" if \code{signif.stars} is
  \code{TRUE}.

  When the \code{"summary.gnm"} object has a \code{"correlation"}
  component, the lower triangle of this matrix is also printed, to two
  decimal places (or symbolically); to see the full matrix of
  correlations print \code{summary(object, correlation =
    TRUE)$correlation} directly.
  
  The standard errors returned by \code{summary.gnm} are scaled by
  \code{sqrt(dispersion)}. If the dispersion is not specified, it is
  taken as \code{1} for the \code{binomial} and \code{Poisson} families,
  and otherwise estimated by the residual Chi-squared statistic divided
  by the residual degrees of freedom. For coefficients that have been
  constrained or are not estimable, the standard error is returned as
  \code{NA}.    
}
\value{
  \code{summary.gnm} returns an object of class \code{"summary.gnm"},
  which is a list with components
  \item{call }{ the \code{"call"} component from object.  }
  \item{ofInterest }{ the \code{"ofInterest"} component from object. }
  \item{family }{ the \code{"family"} component from object. }
  \item{deviance }{ the \code{"deviance"} component from object. }
  \item{aic }{ the \code{"aic"} component from object. }
  \item{df.residual }{ the \code{"df.residual"} component from object. }
  \item{iter }{ the \code{"iter"} component from object. }
  \item{deviance.resid }{ the deviance residuals, see
    \code{\link{residuals.glm}}. } 
  \item{coefficients }{ the matrix of coefficients, standard errors,
    z-values and p-values. }
  \item{elim.coefs }{ if \code{with.eliminate = TRUE} a matrix of
    eliminated coefficients, standard errors, z-values and p-values. }
  \item{dispersion }{ either the supplied argument or the estimated dispersion if
          the latter is \code{NULL}. }
  \item{df}{ a 3-vector of the rank of the model; the number of
    residual degrees of freedom, and number of unconstrained
    coefficients. }
  \item{cov.scaled }{ the estimated covariance matrix scaled by
    \code{dispersion} (see \code{\link{vcov.gnm}} for more details). }
  \item{correlation }{ (only if \code{correlation} is \code{TRUE}) the
    estimated correlations of the estimated coefficients. }
  \item{symbolic.cor }{ (only if \code{correlation} is \code{TRUE}) the value
    of the argument \code{symbolic.cor}. }
}
\note{ The \code{gnm} class includes generalized linear models, and it
  should be noted that \code{summary.gnm} differs from
  \code{\link{summary.glm}} in that it does not omit coefficients which
  are \code{NA} from the objects it returns. (Such coefficients are
  \code{NA} since they have been fixed at \code{0} either by use of the
  \code{constrain} argument to \code{gnm} or by a convention to handle
  linear aliasing).
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
