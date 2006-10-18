\name{nonlin.function}
\alias{nonlin.function}
\title{ Functions to Specify Nonlinear Terms in gnm Models  }
\description{
  Nonlinear terms maybe be specified in the formula argument to gnm by
  a call to a function of class \code{"nonlin"}.  A \code{"nonlin"}
  function takes a list of arguments and returns a list of information
  on the term. 
}
\arguments{
  \item{...}{ arguments required to define the term, e.g. symbolic
    representations of predictors in the term. }
  \item{inst}{(optional) an integer specifying the instance number of
    the term - for compatibility with \code{\link{instances}}. }
}
\value{
  The function should return a list with the following components:
  \item{predictors}{ a list of symbolic expressions or formulae with no
    left hand side which represent (possibly nonlinear)
    predictors that form part of the term. Intercepts will be added by
    default to predictors specified by formulae. If predictors are
    named, these names will be used as a prefix for parameter labels or
    the parameter label itself in the single parameter case (in either
    case, prefixed by the call if supplied.) Matched predictors should
    be given the same name as the matching argument of \code{call}. }
  \item{variables}{ an optional list of expressions representing variables
    in the term. }
  \item{term}{ a function which takes the arguments \code{predLabels}
    and \code{varLabels}, which are vectors of labels defined by
    \code{gnm} that correspond to the specified predictors and variables,
    and returns a deparsed mathematical expression of the full
    term. Only functions recognised by \code{deriv} should be used in
    the expression, e.g. \code{+} rather than \code{sum}.}  
  \item{common}{ an optional numeric index of \code{predictors} with
    duplicated indices identifying single factor predictors for which
    homologous effects are to be estimated. }
  \item{call}{ an optional call to be used as a prefix for parameter labels. }
  \item{match}{ (if \code{call} is non-\code{NULL}) a numeric index
    of \code{predictors} specifying which arguments of \code{call} the
    predictors match to - zero indicating no match. If \code{NULL},
    predictors will be matched sequentially to the arguments of
    \code{call}. It is recommended that matches are specified
    wherever possible, to ensure parameter labels are
    well-defined. Parameters in matched predictors 
    are labelled using "dot-style" labelling, see examples.} 
  \item{start}{ an optional function which takes a named vector of parameters
    corresponding to the predictors and returns a vector of
    starting values for those parameters. This function is ignored if the
    term is nested within another nonlinear term.}
  }
\author{ Heather Turner }
\seealso{ \code{\link{Const}} to specify a constant,
  \code{\link{Dref}} to specify a diagonal reference term,
  \code{\link{Exp}} to specify the exponential of a predictor,
  \code{\link{Inv}} to specify the reciprocal of a predictor,
  \code{\link{Log}} to specify the natural logarithm of a predictor,
  \code{\link{Logit}} to specify the logit of a predictor,
  \code{\link{Mult}} to specify a multiplicative interaction,
  \code{\link{MultHomog}} to specify a homogeneous multiplicative
  interaction,
  \code{\link{Raise}} to raise a predictor to a constant power. }
\examples{
### Equivalent of weighted.MM function in ?nls
weighted.MM <- function(resp, conc){
    list(predictors = list(Vm = substitute(conc), K = 1),
         variables = list(substitute(resp), substitute(conc)),
         term = function(predictors, variables) {
             pred <- paste("(", predictors[1], "/(", predictors[2],
                           " + ", variables[2], "))", sep = "")
             pred <- paste("(", variables[1], " - ", pred, ")/sqrt(",
                           pred, ")", sep = "")
         })
}
class(weighted.MM) <- "nonlin"

## use to fitted weighted Michaelis-Menten model
Treated <- Puromycin[Puromycin$state == "treated", ]
Pur.wt.2 <- gnm( ~ -1 + weighted.MM(rate, conc), data = Treated,
                start = c(Vm = 200, K = 0.1), verbose = FALSE)
Pur.wt.2
## 
## Call:
## gnm(formula = ~-1 + weighted.MM(rate, conc), data = Treated, 
##     start = c(Vm = 200, K = 0.1), verbose = FALSE)
## 
## 
## Coefficients:
##        Vm          K  
## 206.83477    0.05461  
## 
## Deviance:            14.59690 
## Pearson chi-squared: 14.59690 
## Residual df:         10

### The definition of MultHomog
data(occupationalStatus)
MultHomog <- function(..., inst = NULL){
    dots <- match.call(expand.dots = FALSE)[["..."]]
    list(predictors = dots,
         common = rep(TRUE, length(dots)),
         term = function(predictors, ...) {
             paste("(", paste(predictors, collapse = ")*("), ")", sep = "")
         },
         call = as.expression(match.call()))
}
class(MultHomog) <- "nonlin"
## use to fit homogeneous multiplicative interaction
set.seed(1)
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), ofInterest = "[|]",
               family = poisson, data = occupationalStatus,
               verbose = FALSE)
RChomog
## 
## Call:
## 
## gnm(formula = Freq ~ origin + destination + Diag(origin, destination) + 
##     MultHomog(origin, destination), ofInterest = "[|]", family = poisson, 
##     data = occupationalStatus, verbose = FALSE)
## 
## 
## Coefficients of interest:
## MultHomog(., .).`origin|destination`1  
##                              -1.50089  
## MultHomog(., .).`origin|destination`2  
##                              -1.28260  
## MultHomog(., .).`origin|destination`3  
##                              -0.68443  
## MultHomog(., .).`origin|destination`4  
##                              -0.10055  
## MultHomog(., .).`origin|destination`5  
##                              -0.08338  
## MultHomog(., .).`origin|destination`6  
##                               0.42838  
## MultHomog(., .).`origin|destination`7  
##                               0.84452  
## MultHomog(., .).`origin|destination`8  
##                               1.08809  
## 
## Deviance:            32.56098 
## Pearson chi-squared: 31.20716 
## Residual df:         34 
##

## the definition of Exp
Exp <- function(expression, inst = NULL){
    list(predictors = list(substitute(expression)),
         term = function(predictors, ...) {
             paste("exp(", predictors, ")", sep = "")
         },
         call = as.expression(match.call()))
}
class(Exp) <- "nonlin"


## use to fit exponentional model
x <- 1:100
y <- exp(- x / 10)
set.seed(4)
exp1 <- gnm(y ~ Exp(1 + x), verbose = FALSE)
exp1
## 
## Call:
## gnm(formula = y ~ Exp(1 + x), verbose = FALSE)
## 
## 
## Coefficients:
##            (Intercept)  Exp(. + x).(Intercept)            Exp(1 + .).x  
##              1.549e-11              -7.934e-11              -1.000e-01  
## 
## Deviance:            9.342418e-20 
## Pearson chi-squared: 9.342418e-20 
## Residual df:         97
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
