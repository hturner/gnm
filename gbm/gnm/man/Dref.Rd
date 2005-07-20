\name{Dref}
\alias{Dref}
\title{gnm Plug-in Function to Fit Diagonal Reference Terms}
\description{
  A function to provide the objects and functions required to fit a diagonal
  reference term in a generalized nonlinear model using \code{\link{gnm}}.  
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
  \deqn{w_f\gamma_l}{w_f gamma_l}
  where \eqn{w_f} is the weight for factor \eqn{f},
  \eqn{\gamma_l}{gamma_l} is the diagonal effect for level \eqn{l} and
  \eqn{l} is the level of factor \eqn{f} for the given data point.

  The weights are constrained to be nonnegative and to sum to one as
  follows
  \deqn{w_f = \frac{e^{\delta_f}}{\sum_i e^{\delta_i}}}{
    w_f = exp(delta_f)/sum_i(exp(delta_i))}
  and the \eqn{\delta_f}{delta_f} are modelled as specified by the
  \code{formula} argument (constant weights by default). The
  returned parameters are those in the model for
  \eqn{\delta_f}{delta_f}, rather than the implied weights \eqn{w_f}.
  
  To specify a diagonal reference term in the formula argument to \code{gnm},
  the symbolic wrapper \code{Nonlin} must be used, with a call to \code{Dref}
  as the first argument.

  If the factors passed to \code{Dref} do not have exactly the same
  levels, the set of levels in the diagonal reference term is taken to
  be the union of the factor levels, sorted into increasing order.
}
\value{
  A list with the anticipated components of a \code{gnm} plug-in
  function:
  \item{start}{ a vector of starting values for the parameters in the
    diagonal reference term. }
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

  Clifford, P. and Heath, A. F. (1993) The Political Consequences of
  Social Mobility. \emph{J. Roy. Stat. Soc. A}, \bold{156(1)}, 51-61.

  Van der Slik, F. W. P., De Graaf, N. D and Gerris, J. R. M. (2002)
  Conformity to Parental Rules: Asymmetric Influences of Father's and
  Mother's Levels of Education. \emph{European Sociological Review}
  \bold{18(4)}, 489 -- 502. 
}
\author{ Heather Turner }
\seealso{\code{\link{gnm}} and \code{\link{Nonlin}} for usage.
  
  \code{\link{MultHomog}} and \code{\link{LogExcess}} for other
  \code{gnm} plug-in functions.
  
  \code{\link{Mult}} for specifying multiplicative interactions
  in \code{gnm} formulae.}
\examples{
### Examples from Clifford and Heath paper
### (Results differ slightly - possible transcription error in
### published data?)
set.seed(1)
data(voting)
## reconstruct counts voting Labour/non-Labour
count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)

## fit diagonal reference model with constant weights
classMobility <- gnm(yvar ~ Nonlin(Dref(origin, destination)), 
                     family = binomial, data = voting)
prop.table(exp(coef(classMobility)[2:3]))

## create factors indicating movement in and out of salariat (class 1)
upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

## fit separate weights for the "socially mobile" groups
socialMobility <- gnm(yvar ~ Nonlin(Dref(origin, destination,
                                         formula = ~ 1 + downward + upward)),
                      family = binomial, data = voting)
prop.table(exp(coef(socialMobility)[c(4, 7)] + coef(socialMobility)[c(2, 5)]))
prop.table(exp(coef(socialMobility)[c(3, 6)] + coef(socialMobility)[c(2, 5)]))
prop.table(exp(coef(socialMobility)[c(2, 5)]))

## fit separate weights for downwardly mobile groups only
downwardMobility <- gnm(yvar ~ Nonlin(Dref(origin, destination,
                                           formula = ~ 1 + downward)),
                        family = binomial, data = voting)
downwardMobility
prop.table(exp(coef(downwardMobility)[c(3, 5)] + 
               coef(downwardMobility)[c(2, 4)]))
	       prop.table(exp(coef(downwardMobility)[c(2, 4)]))

\dontrun{	       
### Examples from Van der Slik paper
### For illustration only - data not publically available
### Using data in data.frame named 'conformity', with variables
### MCFM - mother's conformity score
### FCFF - father's conformity score
### MOPLM - a factor describing the mother's education with 7 levels
### FOPLF - a factor describing the father's education with 7 levels
### AGEM - mother's birth cohort
### MRMM - mother's traditional role model
### FRMF - father's traditional role model
### MWORK - mother's employment
### MFCM - mother's family conflict score
### FFCF - father's family conflict score

set.seed(1)

## Models for mothers' conformity score as specified in Figure 1
A <- gnm(MCFM ~ -1 + AGEM + MRMM + FRMF + MWORK + MFCM + 
         Nonlin(Dref(MOPLM, FOPLF)), family = gaussian, data = conformity,
         verbose = FALSE)
A
## Call:
## gnm(formula = MCFM ~ -1 + AGEM + MRMM + FRMF + MWORK + MFCM + 
##     Nonlin(Dref(MOPLM, FOPLF)), family = gaussian, data = conformity,
##     verbose = FALSE)
##
##
## Coefficients:
##                    AGEM                      MRMM                      FRMF  
##                 0.06364                  -0.32425                  -0.25324  
##                   MWORK                      MFCM  Dref(MOPLM, FOPLF).MOPLM  
##                -0.06430                  -0.06043                   0.34389  
## Dref(MOPLM, FOPLF).FOPLF      Dref(MOPLM, FOPLF).1      Dref(MOPLM, FOPLF).2  
##                 0.65611                   4.95123                   4.86328  
##     Dref(MOPLM, FOPLF).3      Dref(MOPLM, FOPLF).4      Dref(MOPLM, FOPLF).5  
##                 4.86458                   4.72342                   4.43516  
##     Dref(MOPLM, FOPLF).6      Dref(MOPLM, FOPLF).7  
##                 4.18873                   4.43379  
##
## Deviance:            425.3389 
## Pearson chi-squared: 425.3389 
## Residual df:         576

prop.table(exp(coef(A)[6:7])) ## weights as in Table 4
## Dref(MOPLM, FOPLF).MOPLM Dref(MOPLM, FOPLF).FOPLF 
##               0.4225734                0.5774266

F <- gnm(MCFM ~ -1 + AGEM + MRMM + FRMF + MWORK + MFCM + 
         Nonlin(Dref(MOPLM, FOPLF, formula = ~1 + MFCM)), family = gaussian,
         data = conformity, verbose = FALSE)
F	 
## Call:
## gnm(formula = MCFM ~ -1 + AGEM + MRMM + FRMF + MWORK + MFCM + 
##     Nonlin(Dref(MOPLM, FOPLF, formula = ~1 + MFCM)), family = gaussian, 
##     data = conformity, verbose = FALSE)
##
##
## Coefficients:
##                                                      AGEM  
##                                                   0.05818  
##                                                      MRMM  
##                                                  -0.32701  
##                                                      FRMF  
##                                                  -0.25772  
##                                                     MWORK  
##                                                  -0.07847  
##                                                      MFCM  
##                                                  -0.01694  
## Dref(MOPLM, FOPLF, formula = ~1 + MFCM).MOPLM.(Intercept)  
##                                                   1.03516  
##        Dref(MOPLM, FOPLF, formula = ~1 + MFCM).MOPLM.MFCM  
##                                                  -1.77703  
## Dref(MOPLM, FOPLF, formula = ~1 + MFCM).FOPLF.(Intercept)  
##                                                  -0.03516  
##        Dref(MOPLM, FOPLF, formula = ~1 + MFCM).FOPLF.MFCM  
##                                                   2.77703  
##                 Dref(MOPLM, FOPLF, formula = ~1 + MFCM).1  
##                                                   4.82477  
##                 Dref(MOPLM, FOPLF, formula = ~1 + MFCM).2  
##                                                   4.88066  
##                 Dref(MOPLM, FOPLF, formula = ~1 + MFCM).3  
##                                                   4.83969  
##                 Dref(MOPLM, FOPLF, formula = ~1 + MFCM).4  
##                                                   4.74849  
##                 Dref(MOPLM, FOPLF, formula = ~1 + MFCM).5  
##                                                   4.42019  
##                 Dref(MOPLM, FOPLF, formula = ~1 + MFCM).6  
##                                                   4.17956  
##                 Dref(MOPLM, FOPLF, formula = ~1 + MFCM).7  
##                                                   4.40819  
## 
## Deviance:            420.9022 
## Pearson chi-squared: 420.9022 
## Residual df:         575

prop.table(exp(coef(F))[c(6,8)]) ## weights for MFCM = 0
## Dref(MOPLM, FOPLF, formula = ~1 + MFCM).MOPLM.(Intercept) 
##                                                 0.7446585 
## Dref(MOPLM, FOPLF, formula = ~1 + MFCM).FOPLF.(Intercept) 
##                                                 0.2553415 
prop.table(exp(coef(F)[c(7,9)] + coef(F)[c(6,8)])) ## weights for MFCM = 1
## Dref(MOPLM, FOPLF, formula = ~1 + MFCM).MOPLM.MFCM 
##                                         0.02977851 
## Dref(MOPLM, FOPLF, formula = ~1 + MFCM).FOPLF.MFCM 
##                                         0.97022149 
}
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
