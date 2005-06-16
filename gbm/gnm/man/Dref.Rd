\name{Dref}
\alias{Dref}
\title{``Plug-in'' Function to Fit Diagonal Reference Terms}
\description{
  A function to provide the objects and functions required to fit a diagonal
  reference term in a generalized (non)-linear model using \code{\link{gnm}}.  
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
  Sobel (1981, 1985).  Such terms comprise an additive component for
  each factor of the form
  \deqn{w_i\gamma_l}{w_i gamma_l}
  where \eqn{w_i} is the weight for factor \eqn{i},
  \eqn{\gamma_l}{gamma_l} is the diagonal effect for level \eqn{l} and
  \eqn{l} is the level of factor \eqn{i} for the given data point.

  The weights are constrained to be nonnegative and to sum to one as
  follows
  \deqn{w_i = \frac{e^{\delta_i}}{\sum_r e^{\delta_r}}}{
    w_i = exp(delta_i)/sum_r(exp(delta_r))}
  and the \eqn{\delta_i}{delta_i} are modelled as specified by the
  \code{formula} argument (constant weights by default). The
  returned parameters are those in the model for
  \eqn{\delta_i}{delta_i}, rather than the implied weights \eqn{w_i}
  
  To specify a diagonal reference term in the formula argument to \code{gnm},
  the symbolic wrapper \code{Nonlin} must be used, with a call to \code{Dref}
  as the first argument.

  If the factors passed to \code{Dref} do not have exactly the same
  levels, the set of levels in the diagonal reference term is taken to
  be the union of the factor levels, sorted into increasing order.
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
  \emph{American Sociological Review} \bold{50}, 699--712.

  Van der Slik, F. W. P., De Graaf, N. D and Gerris, J. R. M. (2002)
  Conformity to Parental Rules: Asymmetric Influences of Father's and
  Mother's Levels of Education. \emph{European Sociological Review}
  \bold{18(4)}, 489 -- 502. 
}
\author{ Heather Turner }
\seealso{\code{\link{gnm}} and \code{\link{Nonlin}} for usage.
  
  \code{\link{MultHomog}} for another \code{gnm} plug-in function.
  
  \code{\link{Mult}} for specifying multiplicative interactions
  in \code{gnm} formulae.}
\examples{
### Examples from Van der Slik paper
data(conformity)
set.seed(1)

## Models for mothers' conformity score
A <- gnm(MCFM ~ Nonlin(Dref(MOPLM, FOPLF)) + MFCM + AGEM + FRMF +
         MRMM + MWORK - 1, family = gaussian, data = conformity, trace = TRUE)
prop.table(exp(coef(A)[6:7])) # weights as in Table 4

## equivalent models for fathers
A <- gnm(FCFF ~ Nonlin(Dref(FOPLF, MOPLM)) + FFCF + AGEM + FRMF +
         MRMM + MWORK - 1, family = gaussian, data = conformity)
prop.table(exp(coef(A)[6:7])) # weights as in Table 4

F <- gnm(MCFM ~ Nonlin(Dref(MOPLM, FOPLF, formula = ~ 1 + MFCM)) + MFCM +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = conformity)
prop.table(exp(coef(F))[c(6,8)]) # weights for MFCM = 0
prop.table(exp(coef(F)[c(7,9)] + coef(F)[c(6,8)])) # weights for MFCM = 1
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
