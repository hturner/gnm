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
gnm(formula, constrain = NULL, family = gaussian, data = NULL, subset, weights,
    na.action, method = "gnm.fit", offset, start = NULL,
    control = gnm.control(...), model = TRUE, x = FALSE, vcov = FALSE, ...)
}
\arguments{
  \item{formula}{ a symbolic description of the nonlinear predictor. }
  \item{constrain}{ either a list of indices, or a logical vector,
    indicating elements of the parameter vector to set to zero.}
  \item{family}{ a specification of the error distribution and link function
    to be used in the model. This can be a character string naming
    a family function; a family function, or the result of a call
    to a family function. See \code{\link{family}},
    \code{\link{wedderburn}} and \code{\link{multinomial}} for possibilities.
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
  \item{method}{ the method to be used: either \code{"gnm.fit"} to fit the
    model, \code{"coef"} to return a vector of names for the
    coefficients in the model, or \code{"model.frame"} to return the
    model frame. }
  \item{offset}{ this can be used to specify an a priori known component to
    be added to the predictor during fitting. \code{offset} terms
    can be included in the formula instead or as well, and if both
    are specified their sum is used.}
  \item{start}{ a vector of starting values for the parameters in the
    predictor. }
  \item{control}{ a list of parameters for controlling the fitting process. See
    \code{gnm.control} for details.}
  \item{model}{ logical: if \code{TRUE} the model frame is returned. }
  \item{x}{ logical: if \code{TRUE} the local design matrix from the last
    iteration is returned. }
  \item{vcov}{ logical: if \code{TRUE} the variance-covariance matrix of the
    model parameters is returned. }
  \item{\dots}{ further arguments passed to or from other methods. }
}
\details{
  Models for \code{gnm} are specified by giving a symbolic description
  of the nonlinear predictor, of the form \code{response ~ terms}. The
  \code{response} is typically a numeric vector, see later in this
  section for alternatives. The usual symbolic language may be used to
  specify any linear terms, see \code{formula} for details.

  \code{gnm} has the in-built capability to handle multiplicative
  interactions, which can be specified in the model formula using the
  symbolic wrapper \code{Mult}; e.g. \code{Mult(A, B)} specifies a
  multiplicative interaction between factors \code{A} and
  \code{B}. The family of multiplicative interaction models include
  row-column association models for contingency tables (REF),
  log-multiplicative or UNIDIFF models (REF), and GAMMI models (REF).

  Other nonlinear terms may be incorporated in the model via
  plug-in functions that provide the objects required by \code{gnm} to
  fit the desired term. Such terms are specified in the model formula
  using the symbolic wrapper \code{Nonlin};
  e.g. \code{Nonlin(PlugInFunction(A, B))} specifies a term to be fitted
  by the plug-in function \code{PlugInFunction} involving factors
  \code{A} and \code{B}. The \pkg{gnm} package includes plug-in
  functions for multiplicative interactions with homogeneous effects
  (\code{MultHomog}) and diagonal reference terms (\code{Dref}). Users
  may also define their own plug-in functions, see \code{\link{Nonlin}}
  for details.

  For contingency tables, the data may be provided as an object of class
  \code{"table"} from which the frequencies will be extracted to use
  as the response. In this case, the response should be specified as
  \code{Freq} in the model formula. The \code{"predictors"},
  \code{"fitted.values"}, \code{"residuals"}, \code{"prior.weights"},
  \code{"weights"}, \code{"y"} and \code{"offset"} components of
  the returned \code{gnm} fit will be tables with the same format as the
  data.

  For binomial models, the \code{response} may be specified as a factor
  in which the first level denotes failure and all other levels denote
  success, as a two-column matrix with the columns giving the numbers
  of successes and failures, or as a vector of the proportions of
  successes.

  For multinomial models, the \code{response} should be a factor with
  more than two levels. The model is re-expressed as a poisson model and
  \code{gnm} is re-called. See \code{\link{multinomial}} for an example.

  The \code{gnm} fitting process consists of two stages. In the start-up
  iterations, nonlinear parameters are updated one parameter at a time.
  In the main iterations, all the parameters are updated in each
  iteration. See \code{gnm.control} for more details.

  By default, \code{gnm} uses an over-parameterized representation of
  the model that is being fitted. Only minimal identifiability constraints
  are imposed, so that in general a random parameterization is obtained.
  The parameter estimates are ordered so that those for any linear terms
  appear first.
  
  \code{getContrasts} may be used to obtain estimates of specified
  contrasts, if these contrasts are identifiable. In particular,
  \code{getContrasts} may be used to estimate the contrast between the
  first \code{k - 1} levels of a factor and level \code{k}.

  If appropriate constraints are known in advance, or have been
  determined from a \code{gnm} fit, the model may be (re-)fitted using
  the \code{constrain} argument to specify which parameters should be
  set to zero. \code{update} provides a convenient way of
  re-fitting a \code{gnm} model with new constraints.
}
\value{
  If \code{method = "gnm.fit"}, \code{gnm} returns an object of class
  \code{"gnm"}, which inherits first from \code{"glm"} then \code{"lm"},
  and is a list containing the following components:     
  \item{ call } { the matched call. }
  \item{ formula }{ the formula supplied. }
  \item{ constrain }{ a logical vector, indicating any parameters that were
    constrained to zero in the fitting process. }
  \item{ family }{ the \code{family} object used. }
  \item{ prior.weights }{ the case weights initially supplied. }
  \item{ terms }{ the \code{terms} object used. }
  \item{ na.action }{ the \code{na.action} attribute of the model frame }
  \item{ xlevels }{ a record of the levels of the factors used in fitting. }
  \item{ y }{ the response used. }
  \item{ offset }{ the offset vector used. }
  \item{ control }{ the value of the \code{control} argument used. }
  \item{ coefficients }{ a named vector of coefficients. }
  \item{ predictors }{ the fitted values on the link scale. }
  \item{ fitted.values }{ the fitted mean values, obtained by transforming the
    predictors by the inverse of the link function. }
  \item{ deviance }{ up to a constant, minus twice the maximised
    log-likelihood. Where sensible, the constant is chosen so
    that a saturated model has deviance zero. }
  \item{ aic }{ Akaike's \emph{An Information Criterion}, minus twice the
    maximized log-likelihood plus twice the number of coefficients (so assuming
    that the dispersion is known).}
  \item{ iter }{ the number of main iterations.}
  \item{ conv }{ logical indicating whether the main iterations converged. }
  \item{ weights }{ the \emph{working} weights, that is, the weights used in
    the last iteration.}
  \item{ residuals }{ the \emph{working} residuals, that is, the residuals
    from the last iteration. }
  \item{ df.residual }{ the residual degrees of freedom. }
  \item{ rank }{ the numeric rank of the fitted model. }
  \item{ auxiliary }{ a logical vector indicating whether parameters are
    auxiliary parameters in the re-expression of a multinomial model. }
  
  The list may also contain the components \code{model}, \code{x}, or
  \code{vcov}, if requested in the arguments to \code{gnm}.

  If a binomial \code{gnm} model is specified by giving a two-column
  response, the weights returned by \code{prior.weights} are the total
  numbers of cases (factored by the supplied case weights) and the
  component \code{y} of the result is the proportion of successes.

  For multinomial models, the \code{gnm} object will also contain an
  \code{original.call} component with the original call to
  \code{gnm}. All the other components will relate to the secondary call
  to \code{gnm} which is returned in \code{call}. 
  
  The function \code{summary} may be used to obtain and print a summary
  of the results.

  The generic functions \code{formula}, \code{family}, \code{terms},
  \code{coefficients}, \code{fitted.values}, \code{deviance},
  \code{extractAIC}, \code{weights}, \code{residuals},
  \code{df.residual}, \code{model.frame}, \code{model.matrix} and
  \code{vcov} maybe used to extract components from the object returned
  by \code{gnm} or to construct the relevant objects where necessary.

  Note that the generic functions \code{weights} and \code{residuals} do
  not act as straight-forward accessor functions for \code{gnm} objects,
  but return the prior weights and deviance residuals respectively, as
  for \code{glm} objects.  
}
\references{ Cautres, B, Heath, A F and Firth, D (1998).  Class, religion and vote in
Britain and France.  \emph{La Lettre de la Maison Francaise} \bold{8}. }
\author{ Heather Turner, David Firth }
\seealso{
  \code{\link{formula}} for the symbolic language used to specify
  formulae. 

  \code{Diag} and \code{Symm} for specifying special types of
  interaction.
  
  \code{\link{Mult}}, \code{\link{Nonlin}}, \code{\link{Dref}}
  and \code{\link{MultHomog}} for incorporating nonlinear terms in the
  \code{formula} argument to \code{gnm}.

  \code{\link{residuals.glm}} and the generic functions
  \code{\link{coef}}, \code{\link{fitted}}, etc. for extracting
  components from \code{gnm} objects.

  \code{\link{getContrasts}} to estimate (identifiable) contrasts from a
  \code{gnm} model.
}
\examples{
##  Analysis of a 4-way contingency table

set.seed(1)
data(cautres)
print(cautres)

##  Fit a "double UNIDIFF" model with the religion-vote and class-vote
##  interactions both modulated by nonnegative election-specific
##  multipliers.

doubleUnidiff <- gnm(Freq ~ election:vote + election:class:religion
    + Mult(Exp(election - 1), religion:vote - 1) +
    Mult(Exp(election - 1), class:vote - 1), family = poisson, data = cautres)

##  Examine the multipliers of the class-vote log odds ratios

coefs.of.interest <- grep("Mult2.*election", names(coef(doubleUnidiff)))
coef(doubleUnidiff)[coefs.of.interest]
##  Mult2.Factor1.election1 Mult2.Factor1.election2 
##               -0.5724370               0.1092972 
##  Mult2.Factor1.election3 Mult2.Factor1.election4 
##               -0.1230682              -0.2105843

##  Re-parameterize by setting Mult2.Factor1.election4 to zero

getContrasts(doubleUnidiff, coefs.of.interest)
##                            estimate        se
##  Mult2.Factor1.election1 -0.3618399 0.2534762
##  Mult2.Factor1.election2  0.3198951 0.1320034
##  Mult2.Factor1.election3  0.0875308 0.1446842
##  Mult2.Factor1.election4  0.0000000 0.0000000

##  Same thing but with election 1 as reference category:

getContrasts(doubleUnidiff, rev(coefs.of.interest))
##                            estimate        se
##  Mult2.Factor1.election4  0.3618399 0.2534746
##  Mult2.Factor1.election3  0.4493707 0.2473524
##  Mult2.Factor1.election2  0.6817351 0.2401645
##  Mult2.Factor1.election1  0.0000000 0.0000000

## Re-fit model with Mult2.Factor1.election1 set to zero

doubleUnidiffConstrained <-
    update(doubleUnidiff, constrain = coefs.of.interest[1])

##  Examine the multipliers of the class-vote log odds ratios

coef(doubleUnidiffConstrained)[coefs.of.interest]

## ...as using 'getContrasts' (to 4 d.p.).
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
