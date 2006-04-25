\name{gnm}
\alias{gnm}
\title{ Fitting Generalized Nonlinear Models }
\description{
  \code{gnm} fits generalised nonlinear models using an
  over-parameterised representation. \code{gnm} is able to fit models
  incorporating multiplicative interactions as standard and can fit other
  types of nonlinear effects via \dQuote{plug-in} functions (see details).
}
\usage{
gnm(formula, eliminate = NULL, ofInterest = NULL, constrain = numeric(0),
    constrainTo = numeric(length(constrain)), family = gaussian, 
    data = NULL, subset, weights, na.action,  method = "gnmFit", offset,
    start = NULL, tolerance = 1e-6, iterStart = 2, iterMax = 500,
    trace = FALSE, verbose = TRUE, model = TRUE, x = FALSE,
    termPredictors = FALSE, lsMethod = "qr", ...)  
}
\arguments{
  \item{formula}{ a symbolic description of the nonlinear predictor. }
  \item{eliminate}{ a factor to be included as the first term in the
    model. \code{gnm} will exploit the structure of this factor to
    improve computational efficiency. See details. }
  \item{ofInterest}{ optional coefficients of interest, specified by a
    regular expression, a numeric vector of indices, a character vector of
    names, or "pick" to select from a Tk dialog. If missing, it is
    assumed that all non-\code{eliminate}d coefficients are of interest. }
  \item{constrain}{ coefficients to set to zero, specified by a regular
    expression, a numeric vector of indices, a character vector of 
    names, or "pick" to select from a Tk dialog. }
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
    contain \code{NA}s.  The default is first, any
    \code{na.action} attribute of \code{data}; second, any
    \code{na.action} setting of \code{options}, and third,
    \code{na.fail}.}
  \item{method}{ the method to be used: either \code{"gnmFit"} to fit
    the model using the default maximum likelihood algorithm,
    \code{"coefNames"} to return a character vector of names for the
    coefficients in the model, \code{"model.matrix"} to return the model
    matrix, \code{"model.frame"} to return the model frame, or the name
    of a function providing an alternative fitting algorithm. }
  \item{offset}{ this can be used to specify an a priori known component to
    be added to the predictor during fitting. \code{offset} terms
    can be included in the formula instead or as well, and if both
    are specified their sum is used.}
  \item{start}{ a vector of starting values for the parameters in the
    model; if a starting value is \code{NA}, the default starting value
    will be used. Starting values need not be specified for eliminated
    parameters. } 
  \item{tolerance}{ a positive numeric value specifying the tolerance level for
    convergence. }
  \item{iterStart}{ a positive integer specifying the number of start-up iterations
    to perform. }
  \item{iterMax}{ a positive integer specifying the maximum number of main
    iterations to perform. }
  \item{trace}{ a logical value indicating whether the deviance
    should be printed after each iteration. }
  \item{verbose}{ logical: if \code{TRUE} progress indicators are
    printed as the model is fitted, including a diagnostic error message
    if the algorithm restarts. } 
  \item{model}{ logical: if \code{TRUE} the model frame is returned. }
  \item{x}{ logical: if \code{TRUE} the local design matrix from the last
    iteration is returned. }
  \item{termPredictors}{ logical: if \code{TRUE}, a matrix is returned
    with a column for each term in the model, containing the additive
    contribution of that term to the predictor. }
  \item{lsMethod}{ character: must be one of "chol" or "qr".}
  \item{\dots}{ further arguments passed to fitting function. }
}
\details{
  Models for \code{gnm} are specified by giving a symbolic description
  of the nonlinear predictor, of the form \code{response ~ terms}. The
  \code{response} is typically a numeric vector, see later in this
  section for alternatives. The usual symbolic language may be used to
  specify any linear terms, see \code{\link{formula}} for details.

  \code{gnm} has the in-built capability to handle multiplicative
  interactions, which can be specified in the model formula using the
  symbolic wrapper \code{Mult}; e.g. \code{Mult(A, B)} specifies a
  multiplicative interaction between factors \code{A} and
  \code{B}. The family of multiplicative interaction models include
  row-column association models for contingency tables (e.g., Agresti,
  2002, Sec 9.6), log-multiplicative or UNIDIFF models (Erikson and
  Goldthorpe, 1992; Xie, 1992), and GAMMI models (van Eeuwijk, 1995).

  Other nonlinear terms may be incorporated in the model via
  plug-in functions that provide the objects required by \code{gnm} to
  fit the desired term. Such terms are specified in the model formula
  using the symbolic wrapper \code{Nonlin};
  e.g. \code{Nonlin(PlugInFunction(A, B))} specifies a term to be fitted
  by the plug-in function \code{PlugInFunction} involving factors
  \code{A} and \code{B}. The \pkg{gnm} package includes plug-in
  functions for multiplicative interactions with homogeneous effects
  (\code{MultHomog}) and  diagonal reference terms (\code{Dref}). Users
  may also define their own plug-in functions, see \code{\link{Nonlin}}
  for details. 

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
  the number or iterations reaches \code{iterMax}.    The
  \code{lsMethod} argument specifies
  what numerical method is to be used to solve the
  (typically rank-deficient) least squares problem at the heart of the
  \code{gnm} fitting algorithm: the options are 
  direct solution using a QR decomposition (\code{"qr"}), and matrix
  inversion via Cholesky decomposition (\code{"chol"}).  In both cases,
  the design matrix is standardized and regularized (in the
  Levenberg-Marquardt sense) prior to solving.  If \code{lsMethod} is
  left unspecified, the default is \code{"qr"}, unless \code{eliminate}
  is used in which case the default \code{lsMethod} used is \code{"chol"}. 

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
  
  \code{\link{getContrasts}} may be used to obtain estimates of specified
  contrasts, if these contrasts are identifiable. In particular,
  \code{getContrasts} may be used to estimate the contrasts between the
  first level of a factor and the rest, and obtain standard errors.

  If appropriate constraints are known in advance, or have been
  determined from a \code{gnm} fit, the model may be (re-)fitted using
  the \code{constrain} argument to specify coefficients which should be
  set to zero. Constraints should only be specified for non-eliminated
  parameters. \code{\link{update}} provides a convenient way of re-fitting a 
  \code{gnm} model with new constraints.
}
\value{
  If \code{method = "gnmFit"}, \code{gnm} returns \code{NULL} if the
  algorithm has failed and an object of class \code{"gnm"} otherwise. A
  \code{"gnm"} object inherits first from \code{"glm"} then \code{"lm"}
  and is a list containing the following components:     
  \item{ call }{ the matched call. }
  \item{ formula }{ the formula supplied. }
  \item{constrain}{ a numeric vector specifying any coefficients that
    were constrained in the fitting process. }
  \item{constrainTo}{ a numeric vector of the same length as
    \code{constrain} specifying the values which constrained parameters
    were set to. }  
  \item{ family }{ the \code{family} object used. }
  \item{ prior.weights }{ the case weights initially supplied. }
  \item{ terms }{ the \code{terms} object used. }
  \item{ na.action }{ the \code{na.action} attribute of the model frame }
  \item{ xlevels }{ a record of the levels of the factors used in fitting. }
  \item{ y }{ the response used. }
  \item{ offset }{ the offset vector used. }
  \item{ coefficients }{ a named vector of coefficients. }
  \item{ eliminate }{ the number of eliminated parameters. }
  \item{ ofInterest }{ a named numeric vector of indices, or
    \code{NULL}. }
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
  Agresti, A (2002).  \emph{Categorical Data Analysis} (2nd ed.)  New
  York: Wiley.

  Cautres, B, Heath, A F and Firth, D (1998).  Class,
  religion and vote in Britain and France.  \emph{La Lettre de la Maison
    Francaise} \bold{8}.

  Erikson, R and Goldthorpe, J H (1992).  \emph{The Constant Flux}.
  Oxford: Clarendon Press.

  van Eeuwijk, F A (1995).  Multiplicative interaction in generalized
  linear models.  \emph{Biometrics} \bold{51}, 1017-1032.

  Xie, Y (1992).  The log-multiplicative layer effect model for comparing
  mobility tables.  \emph{American Sociological Review} \bold{57}, 380-395.
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
  
  \code{\link{Mult}}, \code{\link{Nonlin}}, \code{\link{Dref}} and 
  \code{\link{MultHomog}} for incorporating nonlinear terms in the
  \code{formula} argument to \code{gnm}.
  
  \code{\link{residuals.glm}} and the generic functions
  \code{\link{coef}}, \code{\link{fitted}}, etc. for extracting
  components from \code{gnm} objects.

  \code{\link{exitInfo}} to print more information on last iteration
  when \code{gnm} has not converged.

  \code{\link{getContrasts}} to estimate (identifiable) contrasts from a
  \code{gnm} model.
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
                     + Mult(Exp(election - 1), religion:vote - 1) +
                     Mult(Exp(election - 1), class:vote - 1), family = poisson,
                     data = cautres)

##  Examine the multipliers of the class-vote log odds ratios
coefs.of.interest <- grep("Mult2.*election", names(coef(doubleUnidiff)))
coef(doubleUnidiff)[coefs.of.interest]
##  Mult2.Factor1.election1 Mult2.Factor1.election2 
##               -0.5724370               0.1092972 
##  Mult2.Factor1.election3 Mult2.Factor1.election4 
##               -0.1230682              -0.2105843

##  Re-parameterize by setting Mult2.Factor1.election1 to zero
getContrasts(doubleUnidiff, coefs.of.interest)
## [[1]]
##                          estimate        SE    quasiSE    quasiVar
## Mult2.Factor1.election1 0.0000000 0.0000000 0.22854380 0.052232270
## Mult2.Factor1.election2 0.6817370 0.2401642 0.07395880 0.005469905
## Mult2.Factor1.election3 0.4493740 0.2473519 0.09475932 0.008979329
## Mult2.Factor1.election4 0.3618262 0.2534753 0.10934823 0.011957035

##  Same thing but with election 4 as reference category:
getContrasts(doubleUnidiff, rev(coefs.of.interest))
## [[1]]
##                            estimate        SE    quasiSE    quasiVar
## Mult2.Factor1.election4  0.00000000 0.0000000 0.10934823 0.011957035
## Mult2.Factor1.election3  0.08754785 0.1446834 0.09475932 0.008979329
## Mult2.Factor1.election2  0.31991082 0.1320023 0.07395880 0.005469905
## Mult2.Factor1.election1 -0.36182617 0.2534753 0.22854380 0.052232270

##  Re-fit model with Mult2.Factor1.election1 set to zero
doubleUnidiffConstrained <-
    update(doubleUnidiff, constrain = coefs.of.interest[1])

##  Examine the multipliers of the class-vote log odds ratios
coef(doubleUnidiffConstrained)[coefs.of.interest]
##  ...as using 'getContrasts' (to 4 d.p.).
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
