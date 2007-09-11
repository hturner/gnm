\name{Dref}
\alias{Dref}
\alias{DrefWeights}
\title{Specify a Diagonal Reference Term in a gnm Model Formula}
\description{
  Dref is a function of class \code{"nonlin"} to specify a diagonal
  reference term in the formula argument to \code{\link{gnm}}.
}
\usage{
Dref(..., delta = ~ 1)
}
\arguments{
  \item{\dots}{a comma-separated list of two or more factors.}
  \item{delta}{a formula with no left-hand-side specifying the model for
    each factor weight.} 
}
\details{
  \code{Dref} specifies diagonal reference terms as introduced by
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
  \code{delta} argument (constant weights by default). The
  returned parameters are those in the model for
  \eqn{\delta_f}{delta_f}, rather than the implied weights
  \eqn{w_f}. The \code{DrefWeights} function will take a fitted gnm
  model and return the weights \eqn{w_f}, along with their standard
  errors. 

  If the factors passed to \code{Dref} do not have exactly the same
  levels, the set of levels in the diagonal reference term is taken to
  be the union of the factor levels, sorted into increasing order.
}
\value{
  A list with the anticipated components of a "nonlin" function:
  \item{ predictors }{ the factors passed to \code{Dref} and the
    formulae for the weights. }
  \item{ common }{ an index to specify that common effects are to be
    estimated across the factors. }
  \item{ term }{ a function to create a deparsed mathematical expression
    of the term, given labels for the predictors.}
  \item{start }{ a function to generate starting values for the parameters.}
  \item{ call }{ the call to use as a prefix for parameter labels. }
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
\seealso{\code{\link{gnm}}, \code{\link{formula}},
  \code{\link{nonlin.function}}}
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
classMobility <- gnm(yvar ~ -1 + Dref(origin, destination), 
                     family = binomial, data = voting)
DrefWeights(classMobility)

## create factors indicating movement in and out of salariat (class 1)
upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

## fit separate weights for the "socially mobile" groups
socialMobility <- gnm(yvar ~ -1 + Dref(origin, destination,
                                       delta = ~ 1 + downward + upward),
                      family = binomial, data = voting)
DrefWeights(socialMobility)

## fit separate weights for downwardly mobile groups only
downwardMobility <- gnm(yvar ~ -1 + Dref(origin, destination,
                                         delta = ~ 1 + downward),
                        family = binomial, data = voting)
DrefWeights(downwardMobility)

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
         Dref(MOPLM, FOPLF), family = gaussian, data = conformity,
         verbose = FALSE)
A
## Call:
## gnm(formula = MCFM ~ -1 + AGEM + MRMM + FRMF + MWORK + MFCM + 
##     Dref(MOPLM, FOPLF), family = gaussian, data = conformity, 
##     verbose = FALSE)
## 
## Coefficients:
##                     AGEM                      MRMM  
##                  0.06363                  -0.32425  
##                     FRMF                     MWORK  
##                 -0.25324                  -0.06430  
##                     MFCM  Dref(MOPLM, FOPLF)delta1  
##                 -0.06043                  -0.33731  
## Dref(MOPLM, FOPLF)delta2   Dref(., .).MOPLM|FOPLF1  
##                 -0.02505                   4.95121  
##  Dref(., .).MOPLM|FOPLF2   Dref(., .).MOPLM|FOPLF3  
##                  4.86329                   4.86458  
##  Dref(., .).MOPLM|FOPLF4   Dref(., .).MOPLM|FOPLF5  
##                  4.72343                   4.43516  
##  Dref(., .).MOPLM|FOPLF6   Dref(., .).MOPLM|FOPLF7  
##                  4.18873                   4.43378  
## 
## Deviance:            425.3389 
## Pearson chi-squared: 425.3389 
## Residual df:         576    
 
## Weights as in Table 4
DrefWeights(A)
## Refitting with parameters of first Dref weight constrained to zero
## $MOPLM
##    weight        se
## 0.4225636 0.1439829
## 
## $FOPLF
##    weight        se
## 0.5774364 0.1439829 

F <- gnm(MCFM ~ -1 + AGEM + MRMM + FRMF + MWORK + MFCM + 
         Dref(MOPLM, FOPLF, delta = ~1 + MFCM), family = gaussian,
         data = conformity, verbose = FALSE)
F	 
## Call:
## gnm(formula = MCFM ~ -1 + AGEM + MRMM + FRMF + MWORK + MFCM + 
##     Dref(MOPLM, FOPLF, delta = ~1 + MFCM), family = gaussian, 
##     data = conformity, verbose = FALSE)
## 
## 
## Coefficients:
##                                                     AGEM  
##                                                  0.05818  
##                                                     MRMM  
##                                                 -0.32701  
##                                                     FRMF  
##                                                 -0.25772  
##                                                    MWORK  
##                                                 -0.07847  
##                                                     MFCM  
##                                                 -0.01694  
## Dref(MOPLM, FOPLF, delta = ~ . + MFCM).delta1(Intercept)  
##                                                  1.03515  
##           Dref(MOPLM, FOPLF, delta = ~ 1 + .).delta1MFCM  
##                                                 -1.77756  
## Dref(MOPLM, FOPLF, delta = ~ . + MFCM).delta2(Intercept)  
##                                                 -0.03515  
##           Dref(MOPLM, FOPLF, delta = ~ 1 + .).delta2MFCM  
##                                                  2.77756  
##              Dref(., ., delta = ~ 1 + MFCM).MOPLM|FOPLF1  
##                                                  4.82476  
##              Dref(., ., delta = ~ 1 + MFCM).MOPLM|FOPLF2  
##                                                  4.88066  
##              Dref(., ., delta = ~ 1 + MFCM).MOPLM|FOPLF3  
##                                                  4.83969  
##              Dref(., ., delta = ~ 1 + MFCM).MOPLM|FOPLF4  
##                                                  4.74850  
##              Dref(., ., delta = ~ 1 + MFCM).MOPLM|FOPLF5  
##                                                  4.42020  
##              Dref(., ., delta = ~ 1 + MFCM).MOPLM|FOPLF6  
##                                                  4.17957  
##              Dref(., ., delta = ~ 1 + MFCM).MOPLM|FOPLF7  
##                                                  4.40819  
## 
## Deviance:            420.9022 
## Pearson chi-squared: 420.9022 
## Residual df:         575 
##
##

## Standard error for MFCM == 1 lower than reported by Van der Slik et al
DrefWeights(F)
## Refitting with parameters of first Dref weight constrained to zero
## $MOPLM
##   MFCM     weight        se
## 1    1 0.02974675 0.2277711
## 2    0 0.74465224 0.2006916
## 
## $FOPLF
##   MFCM    weight        se
## 1    1 0.9702532 0.2277711
## 2    0 0.2553478 0.2006916
}
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
