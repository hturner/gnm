\name{LogExcess}
\alias{LogExcess}
\title{ gnm Plug-in Function to Fit Log Excess Terms  }
\description{
  A function to provide the objects and functions required to fit log
  excess terms in a generalized nonlinear model using \code{\link{gnm}}.
}
\usage{
LogExcess(x, side = "left",
          constraint = ifelse(side == "right", max(x) + 1, min(x) - 1))
}
\arguments{
  \item{x}{ a numeric vector }
  \item{side}{ either "left" or "right" to indicate the side on which
    the excess is taken }
  \item{constraint}{ a constraint on the excess parameter: for \code{side ==
    "left"}, this can be at most \code{min(x)}; for \code{side ==
    "right"}, at least \code{max(x)}. }
}
\details{
  For \code{side == "left"}, the \code{LogExcess} term is defined as
  \deqn{\beta_1\log(x - a) = \beta_1\log(x - (u - \exp(a*))}{
   beta_1 log(x - a) = beta_1 log(x - (u - exp(a*)))}
  where the excess parameter \code{a} is constrained by the upper limit
  \code{u}. \code{gnm} estimates the parameters \eqn{\beta_1}{beta_1}
  and \eqn{a*}.

  For \code{side == "right"}, the \code{LogExcess} term is defined as
  \deqn{\beta_2\log(b - x) = \beta_2\log((l + \exp(b*)) - x}{
   beta_2 log(b - x) = beta_2 log((l + exp(b*)) - x)}
  where the excess parameter \code{b} is constrained by the lower limit
  \code{l}.\code{gnm} estimates the parameters \eqn{\beta_2}{beta_2}
  and \eqn{b*}.
}
\value{
  A list with the anticipated components of a \code{gnm} plug-in
  function:
  \item{starting values}{ a vector of starting values for the parameters
  in the log excess term }
  \item{labels}{ a character vector of labels for the parameters in the
    log excess term. }
  \item{predictor}{ a function that takes estimates of the parameters in the
    log excess term and returns the fitted values. }
  \item{localDesignFunction}{ a function that takes a vector of estimates of
    the parameters in the log excess term and returns
    the local design matrix. }
}
\references{  }
\author{ Heather Turner }
\seealso{ \code{\link{gnm}} and \code{\link{Nonlin}} for usage.
  
  \code{\link{Dref}} and \code{\link{MultHomog}} for other \code{gnm}
  plug-in functions. }
\examples{

}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
