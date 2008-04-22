\name{gnm-package}
\alias{gnm-package}
\docType{package}
\title{
Generalized Nonlinear Models
}
\description{
Functions to specify, fit and evaluate generalized nonlinear
models.
}
\details{
  \code{gnm} provides functions to fit generalized nonlinear models by
  maximum likelihood. Such models extend the class of generalized linear
  models by allowing nonlinear terms in the predictor.

  Some special cases are models with multiplicative interaction terms,
  such as the UNIDIFF and row-column association models from sociology
  and the AMMI and GAMMI models from crop science; stereotype models for
  ordered categorical response, and diagonal reference models for
  dependence on a square two-way classification.
  
  \code{gnm} is a major re-working of an earlier Xlisp-Stat package,
  "Llama". Over-parameterized representations of models are used
  throughout; functions are provided for inference on estimable
  parameter combinations, as well as standard methods for diagnostics
  etc.

  The following documentation provides further information on the
  \code{gnm} package:
  \describe{
    \item{gnmOverview}{\code{vignette("gnmOverview", package = "gnm")}}
    \item{CHANGES}{\code{file.show(system.file("CHANGES", package = "gnm"))}}
  }
}
\author{
Heather Turner and David Firth

Maintainer: Heather Turner <h.turner@warwick.ac.uk>
}
\references{
  http://www.warwick.ac.uk/go/heatherturner/gnm
}
\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\seealso{
  \code{\link{gnm}} for the model fitting function, with links to
  associated functions.
}
\examples{
demo(gnm)
}
