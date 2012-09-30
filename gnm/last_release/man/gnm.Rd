\name{gnm}
\alias{gnm}
\title{ Fitting Generalized Nonlinear Models }
\description{
  \code{gnm} fits generalised nonlinear models using an
  over-parameterised representation. Nonlinear terms are specified by
  calls to functions of class \code{"nonlin"}.
}
\usage{
gnm(formula, eliminate = NULL, ofInterest = NULL, constrain = numeric(0),
    constrainTo = numeric(length(constrain)), family = gaussian, 
    data = NULL, subset, weights, na.action,  method = "gnmFit", 
    checkLinear = TRUE, offset, start = NULL, etastart = NULL,
    mustart = NULL, tolerance = 1e-06, iterStart = 2, iterMax = 500,
    trace = FALSE, verbose = TRUE, model = TRUE, x = TRUE,
    termPredictors = FALSE, ridge = 1e-08, ...) 
}
\arguments{
  \item{formula}{ a symbolic description of the nonlinear predictor. }
  \item{eliminate}{ a factor to be included as the first term in the
    model. \code{gnm} will exploit the structure of this factor to
    improve computational efficiency. See details. }
  \item{ofInterest}{ optional coefficients of interest, specified by a
    regular expression, a numeric vector of indices, a character vector of
    names, or "[?]" to select from a Tk dialog. If missing, it is
    assumed that all non-\code{eliminate}d coefficients are of interest. }
  \item{constrain}{ (non-eliminated) coefficients to constrain, specified
    by a regular expression, a numeric vector of indices, a logical vector, a
    character vector of names, or "[?]" to select from a Tk dialog. }
  \item{constrainTo}{ a numeric vector of the same length as
    \code{constrain} specifying the values to constrain to. By default
    constrained parameters will be set to zero. } 
  \item{family}{ a specification of the error distribution and link function
    to be used in the model. This can be a character string naming
    a family function; a family function, or the result of a call
    to a family function. See \code{\link{family}} and
    \code{\link{wedderburn}} for possibilities. 
  }
  \item{data}{ an optional data frame containing the variables in the model.
    If not found in \code{data}, the variables are taken from
    \code{environment(formula)}, typically the environment from which
    \code{gnm} is called.}
  \item{subset}{ an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{weights}{ an optional vector of weights to be used in the fitting
    process.}
  \item{na.action}{ a function which indicates what should happen when the data
    contain \code{NA}s. If \code{data} is a contingency table, the
    default is \code{"exclude"}. Otherwise the default is first, any
    \code{na.action} attribute of \code{data}; second, any
    \code{na.action} setting of \code{options}, and third,
    \code{na.fail}.}
  \item{method}{ the method to be used: either \code{"gnmFit"} to fit
    the model using the default maximum likelihood algorithm,
    \code{"coefNames"} to return a character vector of names for the
    coefficients in the model, \code{"model.matrix"} to return the model
    matrix, \code{"model.frame"} to return the model frame, or the name
    of a function providing an alternative fitting algorithm. }
  \item{checkLinear}{ logical: if \code{TRUE} \code{glm.fit} is used
    when the predictor is found to be linear }
  \item{offset}{ this can be used to specify an a priori known component to
    be added to the predictor during fitting. \code{offset} terms
    can be included in the formula instead or as well, and if both
    are specified their sum is used.}
  \item{start}{ a vector of starting values for the parameters in the
    model; if a starting value is \code{NA}, the default starting value
    will be used. Starting values need not be specified for eliminated
    parameters. }
  \item{etastart}{ starting values for the linear predictor. }
  \item{mustart}{ starting values for the vector of means. }
  \item{tolerance}{ a positive numeric value specifying the tolerance level for
    convergence. }
  \item{iterStart}{ a positive integer specifying the number of start-up iterations
    to perform. }
  \item{iterMax}{ a positive integer specifying the maximum number of main
    iterations to perform. }
  \item{trace}{ a logical value indicating whether the deviance
    should be printed after each iteration. }
  \item{verbose}{ logical: if \code{TRUE} and model includes nonlinear
    terms, progress indicators are printed as the model is fitted,
    including a diagnostic error message if the algorithm fails. } 
  \item{model}{ logical: if \code{TRUE} the model frame is returned. }
  \item{x}{ logical: if \code{TRUE} the local design matrix from the last
    iteration is included as a component of returned model object. }
  \item{termPredictors}{ logical: if \code{TRUE}, a matrix is returned
    with a column for each term in the model, containing the additive
    contribution of that term to the predictor. }
  \item{ridge}{numeric, a positive value for the ridge constant to be
    used in the fitting algorithm}
  \item{\dots}{ further arguments passed to fitting function. }
}
\details{
  Models for \code{gnm} are specified by giving a symbolic description
  of the nonlinear predictor, of the form \code{response ~ terms}. The
  \code{response} is typically a numeric vector, see later in this
  section for alternatives. The usual symbolic language may be used to
  specify any linear terms, see \code{\link{formula}} for details.

  Nonlinear terms may be specified by calls to functions of class
  "nonlin". There are several "nonlin" functions in the \code{gnm}
  package. Some of these specify simple
  mathematical functions of predictors: \code{Exp}, \code{Mult}, and 
  \code{Inv}. Others specify more specialised nonlinear terms, in
  particular \code{MultHomog} specifies homogeneous multiplicative
  interactions and \code{Dref} specifies diagonal reference terms. Users
  may also define their own "nonlin" functions, see
  \code{\link{nonlin.function}} for details.

  The \code{eliminate} argument may be used to specify a factor that
  is to be included as the first term in the model (since an intercept
  is then redundant, none is fitted). The structure of the factor is
  exploited to improve computational efficiency --- substantially so if
  the \code{eliminate}d factor has a large number of levels. Use of
  \code{eliminate} is designed for factors that are required in the
  model but are not of direct interest (e.g., terms needed to fit
  multinomial-response models as conditional Poisson models). See
  \code{\link{backPain}} for an example.

  The \code{ofInterest} argument may be used to specify coefficients of
  interest, the indices of which are returned in the \code{ofInterest}
  component of the model object. \code{print()} displays of the model
  object or its components obtained using accessor functions such as
  \code{coef()} etc, will only show these coefficients. In addition
  methods for \code{"gnm"} objects which may be applied to a subset of
  the parameters are by default applied to the coefficients of interest.
  See \code{\link{ofInterest}} for accessor and replacement functions.

  For contingency tables, the data may be provided as an object of class
  \code{"table"} from which the frequencies will be extracted to use
  as the response. In this case, the response should be specified as
  \code{Freq} in the model formula. The \code{"predictors"},
  \code{"fitted.values"}, \code{"residuals"}, \code{"prior.weights"},
  \code{"weights"}, \code{"y"} and \code{"offset"} components of
  the returned \code{gnm} fit will be tables with the same format as the
  data, completed with \code{NA}s where necessary.

  For binomial models, the \code{response} may be specified as a factor
  in which the first level denotes failure and all other levels denote
  success, as a two-column matrix with the columns giving the numbers
  of successes and failures, or as a vector of the proportions of
  successes.

  The \code{gnm} fitting algorithm consists of two stages.  In the start-up
  iterations, any nonlinear parameters that are not specified by either the
  \code{start} argument of \code{gnm} or a plug-in function are
  updated one parameter at a time, then the linear parameters are
  jointly updated before the next iteration. In the main iterations, all
  the parameters are jointly updated, until convergence is reached or
  the number or iterations reaches \code{iterMax}. To solve the
  (typically rank-deficient) least squares problem at the heart of the
  \code{gnm} fitting algorithm,  the design matrix is standardized and
  regularized (in the Levenberg-Marquardt sense) prior to solving; the
  \code{ridge} argument provides a degree of control over the
  regularization performed (smaller values may sometimes give faster
  convergence but can lead to numerical instability).

  Convergence is judged by comparing the squared components of the score vector
  with corresponding elements of the diagonal of the Fisher information
  matrix. If, for all components of the score vector, the ratio is less
  than \code{tolerance^2}, or the corresponding diagonal element of the
  Fisher information matrix is less than 1e-20, iterations cease. If the
  algorithm has not converged by \code{iterMax} iterations,
  \code{\link{exitInfo}} can be used to print information on the
  parameters which failed the convergence criteria at the last iteration.

  By default, \code{gnm} uses an over-parameterized representation of
  the model that is being fitted. Only minimal identifiability constraints
  are imposed, so that in general a random parameterization is obtained.
  The parameter estimates are ordered so that those for any linear terms
  appear first.
  
  \code{\link{getContrasts}} may be used to obtain estimates of
  specified scaled contrasts, if these contrasts are identifiable. For
  example, \code{getContrasts} may be used to estimate the contrasts
  between the first level of a factor and the rest, and obtain standard
  errors. 

  If appropriate constraints are known in advance, or have been
  determined from a \code{gnm} fit, the model may be (re-)fitted using
  the \code{constrain} argument to specify coefficients which should be
  set to values specified by \code{constrainTo}. Constraints should only
  be specified for non-eliminated parameters. \code{\link{update}}
  provides a convenient way of re-fitting a  \code{gnm} model with new
  constraints. 
}
\value{
  If \code{method = "gnmFit"}, \code{gnm} returns \code{NULL} if the
  algorithm has failed and an object of class \code{"gnm"} otherwise. A
  \code{"gnm"} object inherits first from \code{"glm"} then \code{"lm"}
  and is a list containing the following components:     
  \item{ call }{ the matched call. }
  \item{ formula }{ the formula supplied. }
  \item{ constrain }{ a numeric vector specifying any coefficients that
    were constrained in the fitting process. }
  \item{ constrainTo }{ a numeric vector of the same length as
    \code{constrain} specifying the values which constrained parameters
    were set to. }  
  \item{ family }{ the \code{family} object used. }
  \item{ prior.weights }{ the case weights initially supplied. }
  \item{ terms }{ the \code{terms} object used. }
  \item{ data }{ the \code{data} argument. }
  \item{ na.action }{ the \code{na.action} attribute of the model frame }
  \item{ xlevels }{ a record of the levels of the factors used in fitting. }
  \item{ y }{ the response used. }
  \item{ offset }{ the offset vector used. }
  \item{ coefficients }{ a named vector of non-eliminated coefficients,
    with an attribute \code{"eliminated"} specifying the eliminated
    coefficients if \code{eliminate} is non-\code{NULL}. }
  \item{ eliminate }{ the \code{eliminate} argument. }
  \item{ ofInterest }{ a named numeric vector of indices corresponding
    to non-eliminated coefficients, or \code{NULL}. }
  \item{ predictors }{ the fitted values on the link scale. }
  \item{ fitted.values }{ the fitted mean values, obtained by transforming the
    predictors by the inverse of the link function. }
  \item{ deviance }{ up to a constant, minus twice the maximised
    log-likelihood. Where sensible, the constant is chosen so
    that a saturated model has deviance zero. }
  \item{ aic }{ Akaike's \emph{An Information Criterion}, minus twice the
    maximized log-likelihood plus twice the number of parameters (so assuming
    that the dispersion is known).}
  \item{ iter }{ the number of main iterations.}
  \item{ conv }{ logical indicating whether the main iterations
    converged, with an attribute for use by \code{\link{exitInfo}} if
    \code{FALSE}. } 
  \item{ weights }{ the \emph{working} weights, that is, the weights used in
    the last iteration.}
  \item{ residuals }{ the \emph{working} residuals, that is, the residuals
    from the last iteration. }
  \item{ df.residual }{ the residual degrees of freedom. }
  \item{ rank }{ the numeric rank of the fitted model. }
  
  The list may also contain the components \code{model}, \code{x},
  or \code{termPredictors} if requested in the arguments to \code{gnm}.

  If a table was passed to \code{data} and the default for
  \code{na.action} was not overridden, the list will also contain a
  \code{table.attr} component, for use by the extractor functions.

  If a binomial \code{gnm} model is specified by giving a two-column
  response, the weights returned by \code{prior.weights} are the total
  numbers of cases (factored by the supplied case weights) and the
  component \code{y} of the result is the proportion of successes.
  
  The function \code{\link{summary.gnm}} may be used to obtain and print
  a summary of the results, whilst \code{\link{plot.gnm}} may be used
  for model diagnostics.

  The generic functions \code{\link{formula}}, \code{\link{family}}, 
  \code{\link{terms}}, \code{\link{coefficients}}, 
  \code{\link{fitted.values}}, \code{\link{deviance}}, 
  \code{\link{extractAIC}}, \code{\link{weights}}, 
  \code{\link{residuals}}, \code{\link{df.residual}}, 
  \code{\link{model.frame}}, \code{\link{model.matrix}}, 
  \code{\link{vcov}} and \code{\link{termPredictors}} maybe used to
  extract components from the object returned by \code{\link{gnm}} or to
  construct the relevant objects where necessary. 

  Note that the generic functions \code{\link{weights}} and
  \code{\link{residuals}} do not act as straight-forward accessor
  functions for \code{gnm} objects, but return the prior weights and
  deviance residuals respectively, as for \code{glm} objects.  
}
\references{
  Cautres, B, Heath, A F and Firth, D (1998).  Class,
  religion and vote in Britain and France.  \emph{La Lettre de la Maison
    Francaise} \bold{8}.
}
\author{ Heather Turner, David Firth }
\note{
  Regular expression matching is performed using \code{grep} with
  default settings.
}
\seealso{
  \code{\link{formula}} for the symbolic language used to specify
  formulae. 

  \code{\link{Diag}} and \code{\link{Symm}} for specifying special types
  of interaction.
  
  \code{Exp}, \code{Mult}, \code{Inv}, \code{\link{MultHomog}},
  \code{\link{Dref}} and \code{\link{nonlin.function}} for incorporating
  nonlinear terms in the  \code{formula} argument to \code{gnm}.
  
  \code{\link{residuals.glm}} and the generic functions
  \code{\link{coef}}, \code{\link{fitted}}, etc. for extracting
  components from \code{gnm} objects.

  \code{\link{exitInfo}} to print more information on last iteration
  when \code{gnm} has not converged.

  \code{\link{getContrasts}} to estimate (identifiable) scaled contrasts
  from a \code{gnm} model.
}
\examples{
###  Analysis of a 4-way contingency table
set.seed(1)
data(cautres)
print(cautres)

##  Fit a "double UNIDIFF" model with the religion-vote and class-vote
##  interactions both modulated by nonnegative election-specific
##  multipliers.
doubleUnidiff <- gnm(Freq ~ election:vote + election:class:religion
                     + Mult(Exp(election), religion:vote) +
                     Mult(Exp(election), class:vote), family = poisson,
                     data = cautres)

##  Examine the multipliers of the class-vote log odds ratios
ofInterest(doubleUnidiff) <- pickCoef(doubleUnidiff, "class:vote[).]")
coef(doubleUnidiff)
## Coefficients of interest:
## Mult(Exp(.), class:vote).election1 
##                        -0.38357138 
## Mult(Exp(.), class:vote).election2 
##                         0.29816599 
## Mult(Exp(.), class:vote).election3 
##                         0.06580307 
## Mult(Exp(.), class:vote).election4 
##                        -0.02174104

##  Re-parameterize by setting first multiplier to zero
getContrasts(doubleUnidiff, ofInterest(doubleUnidiff))
##                                     estimate        SE
## Mult(Exp(.), class:vote).election1 0.0000000 0.0000000
## Mult(Exp(.), class:vote).election2 0.6817374 0.2401644
## Mult(Exp(.), class:vote).election3 0.4493745 0.2473521
## Mult(Exp(.), class:vote).election4 0.3618301 0.2534754
##                                       quasiSE    quasiVar
## Mult(Exp(.), class:vote).election1 0.22854401 0.052232363
## Mult(Exp(.), class:vote).election2 0.07395886 0.005469913
## Mult(Exp(.), class:vote).election3 0.09475938 0.008979340
## Mult(Exp(.), class:vote).election4 0.10934798 0.011956981

##  Same thing but with last multiplier as reference category:
getContrasts(doubleUnidiff, rev(ofInterest(doubleUnidiff)))
##                                       estimate        SE
## Mult(Exp(.), class:vote).election4  0.00000000 0.0000000
## Mult(Exp(.), class:vote).election3  0.08754436 0.1446833
## Mult(Exp(.), class:vote).election2  0.31990727 0.1320022
## Mult(Exp(.), class:vote).election1 -0.36183013 0.2534754
##                                       quasiSE    quasiVar
## Mult(Exp(.), class:vote).election4 0.10934798 0.011956981
## Mult(Exp(.), class:vote).election3 0.09475938 0.008979340
## Mult(Exp(.), class:vote).election2 0.07395886 0.005469913
## Mult(Exp(.), class:vote).election1 0.22854401 0.052232363

##  Re-fit model with first multiplier set to zero
doubleUnidiffConstrained <-
    update(doubleUnidiff, constrain = ofInterest(doubleUnidiff)[1])

##  Examine the multipliers of the class-vote log odds ratios
coef(doubleUnidiffConstrained)[ofInterest(doubleUnidiff)]
##  ...as using 'getContrasts' (to 4 d.p.).
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
