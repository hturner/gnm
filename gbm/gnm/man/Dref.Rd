\name{Dref}
\alias{Dref}
\title{``Plug-in'' Function to Fit Diagonal Reference Terms}
\description{
  A function to provide the objects and functions required to fit a diagonal
  reference term in a generalized (non)-linear model using \code{gnm}.  
}
\usage{
Dref(..., formula = ~ 1)
}
\arguments{
  \item{\dots}{a comma-separated list of two or more factors.}
  \item{formula}{a symbolic description of the model for each factor weight.}
}
\details{
  Designed to work as a plug-in function for \code{gnm}, \code{Dref} produces
  the objects required to fit diagonal reference terms as introduced by
  Sobel (1981, 1985).  
  
  To specify a diagonal reference term in the formula argument to \code{gnm},
  the symbolic wrapper \code{Nonlin} must be used, with a call to \code{Dref}
  as the first argument.  
  
  By default, \code{Dref} will produce the objects required to fit a diagonal
  reference term involving the factors specified, with a constant weight for
  each factor. An alternative model for the weights may be specified via the
  \code{formula} argument.  
}
\value{
  A list with the components required of a \code{gnm} plug-in function:
  \item{labels}{a character vector of labels for the parameters in the
    diagonal reference term.}
  \item{predictor}{a function that takes estimates of the parameters in the
    diagonal reference term and returns a matrix whose columns
    are the fitted additive components of this term.}
  \item{localDesignFunction}{a function that takes the result of the
    predictor function and returns the local design
    matrix.}
}
\references{Sobel, M. E. (1981), Diagonal mobility models: A substantively
  motivated class of designs for the analysis of mobility effects.
  \emph{American Sociological Review} \bold{46}, 893--906.
  
  Sobel, M. E. (1985), Social mobility and fertility revisited: Some
  new models for the analysis of the mobility effects hypothesis.  
  \emph{American Sociological Review} \bold{50}, 699--712.}
\author{ Heather Turner }
\seealso{\code{\link{gnm}} and \code{\link{Nonlin}} for usage.
  
  \code{\link{MultHomog}} for another \code{gnm} plug-in function.
  
  \code{\link{Mult}} for specifying multiplicative interactions
  in \code{gnm} formulae.}
\examples{
    
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
