\name{getContrasts}
\alias{getContrasts}
\title{ Estimated Contrasts and Standard Errors for Parameters in a gnm
 Model } 
\description{
Computes contrasts or scaled contrasts for a set of (non-eliminated)
parameters from a \code{\link{gnm}} model, and computes standard errors
for the estimated contrasts.  Where possible, quasi standard errors are
also computed. 
}
\usage{
getContrasts(model, set = NULL,  ref = "first", scaleRef = "mean",
  scaleWeights = NULL, dispersion = NULL, use.eliminate = TRUE, ...)
}
\arguments{
  \item{model}{ a model object of class \code{"gnm"}.}
  \item{set}{ a vector of indices (numeric) or coefficient names
    (character). If \code{NULL}, a dialog will open for parameter
    selection. }
  \item{ref}{either a single numeric index, or a vector
    of real numbers which sum to 1, or one of the character
    strings \code{"first"}, \code{"last"} or \code{"mean"}.}
  \item{scaleRef}{as for \code{ref}}
  \item{scaleWeights}{either \code{NULL}, a vector of real numbers,
    \code{"unit"} or \code{"setLength"}.}
  \item{dispersion}{either \code{NULL}, or a positive number by which
    the model's variance-covariance matrix should be scaled.}
  \item{\dots}{ arguments to pass to other functions. }
}
\details{
  The indices in \code{set}
  must all be in \code{1:length(coef(object))}.  If
  \code{set = NULL}, a dialog is presented for the selection
  of indices (model coefficients).

  For the set of coefficients selected, contrasts and their standard
  errors are computed. A check is performed first on the estimability of
  all such contrasts.  The specific contrasts to be computed are
  controlled by the choice of \code{ref}: this may be \code{"first"}
  (the default), for contrasts with the first of the selected
  coefficients, or \code{"last"} for contrasts with the last, or
  \code{"mean"} for contrasts with the arithmetic mean of the
  coefficients in the selected set; or it may be an arbitrary vector of
  weights (summing to 1, not necessarily all non-negative) which specify
  a weighted mean against which contrasts are taken; or it may be a
  single index specifying one of the coefficients with which all
  contrasts should be taken.  Thus, for example, \code{ref = 1} is
  equivalent to \code{ref = "first"}, and \code{ref = c(1/3, 1/3, 1/3)}
  is equivalent to \code{ref = "mean"} when there are three coefficients
  in the selected \code{set}.

  The contrasts may be scaled by
  \deqn{\frac{1}{\sqrt{\sum_r v_r * d_r^2}}}{1/sqrt(sum(v * d))}
  where \eqn{d_r} is a contrast of the r'th coefficient in \code{set}
  with the reference level specified by \code{scaleRef} and \eqn{v} is a
  vector of weights (of the same length as \code{set})
  specified by \code{scaleWeights}. If
  \code{scaleWeights} is \code{NULL} (the default), \code{scaleRef}
  is ignored and no scaling is performed. Other options for
  \code{scaleWeights} are \code{"unit"} for weights equal to one and
  \code{"setLength"} for weights equal to the reciprocal of
  \code{length(set)}. If \code{scaleRef}
  is the same as \code{ref}, these options constrain the sum of squared
  contrasts to 1 and \code{length(set)} respectively.

  Quasi-variances (and corresponding quasi standard errors) are
  reported for \bold{unscaled} contrasts where possible. These statistics are
  invariant to the choice of \code{ref}, see Firth (2003) or Firth and
  Menezes (2004) for more details.
}
\value{
  An object of class
  \code{qv} --- see \code{\link[qvcalc]{qvcalc}}.
}
\author{ David Firth }
\seealso{ \code{\link{gnm}}, \code{\link{se}},
  \code{\link{checkEstimable}}, \code{\link[qvcalc]{qvcalc}},
  \code{\link{ofInterest}}}
\references{
  Firth, D (2003).  Overcoming the reference category problem in the
  presentation of statistical models.  \emph{Sociological Methodology}
  \bold{33}, 1--18.

  Firth, D and Menezes, R X de (2004).  Quasi-variances.
  \emph{Biometrika} \bold{91}, 65--80.
  }
  \examples{
### Unscaled contrasts ###
set.seed(1)
data(yaish)

## Fit the "UNIDIFF" mobility model across education levels -- see ?yaish
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(educ), orig:dest),
               ofInterest = "[.]educ", family = poisson,
               data = yaish,  subset = (dest != 7))
## Examine the education multipliers (differences on the log scale):
unidiffContrasts <- getContrasts(unidiff, ofInterest(unidiff))
plot(unidiffContrasts,
  main = "Unidiff multipliers (log scale): intervals based on
           quasi standard errors",
  xlab = "Education level", levelNames = 1:5)


### Scaled contrasts (elliptical contrasts) ###
set.seed(1)
data(mentalHealth)

##  Goodman Row-Column association model fits well (deviance 3.57, df 8)
mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)
RC1model <- gnm(count ~ SES + MHS + Mult(SES, MHS),
                family = poisson, data = mentalHealth)
## Row scores and column scores are both unnormalized in this
## parameterization of the model 

## The scores can be normalized as in Agresti's eqn (9.15):
rowProbs <- with(mentalHealth, tapply(count, SES, sum) / sum(count))
colProbs <- with(mentalHealth, tapply(count, MHS, sum) / sum(count))
mu <- getContrasts(RC1model, pickCoef(RC1model, "[.]SES"),
                   ref = rowProbs, scaleRef = rowProbs,
                   scaleWeights = rowProbs)
nu <- getContrasts(RC1model, pickCoef(RC1model, "[.]MHS"),
                   ref = colProbs, scaleRef = colProbs,
                   scaleWeights = colProbs)
all.equal(sum(mu$qv[,1] * rowProbs), 0)
all.equal(sum(nu$qv[,1] * colProbs), 0)
all.equal(sum(mu$qv[,1]^2 * rowProbs), 1)
all.equal(sum(nu$qv[,1]^2 * colProbs), 1)
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
