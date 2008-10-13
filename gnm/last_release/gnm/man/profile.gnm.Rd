\name{profile.gnm}
\alias{profile.gnm}
\alias{plot.profile.gnm}
\title{ Profile Deviance for Parameters in a Generalized Nonlinear Model }
\description{
  For one or more parameters in a generalized nonlinear model, profile
  the deviance over a range of values about the fitted estimate.
}
\usage{
\method{profile}{gnm}(fitted, which = ofInterest(fitted), alpha = 0.05, maxsteps = 10,
            stepsize = NULL, trace = FALSE, ...) 
}
\arguments{
  \item{fitted}{ an object of class \code{"gnm"}. }
  \item{which}{ (optional) either a numeric vector of indices or a
    character vector of names, specifying the parameters over which the
    deviance is to be profiled. If missing, the deviance is profiled
    over all parameters. }  
  \item{alpha}{ the significance level of the z statistic, indicating
    the range that the profile must cover (see details). }
  \item{maxsteps}{ the maximum number of steps to take either side of
    the fitted parameter. }
  \item{stepsize}{ (optional) a numeric vector of length two, specifying
    the size of steps to take when profiling down and up respectively,
    or a single number specifying the step size in both directions. If
    missing, the step sizes will be determined automatically. }
  \item{trace}{ logical, indicating whether profiling should be
    traced. } 
  \item{\dots}{ further arguments. }
}
\details{
  This is a method for the generic function \code{profile} in the
  \code{base} package.

  For a given parameter, the deviance is profiled by constraining that
  parameter to certain values either side of its estimate in the fitted
  model and refitting the model.

  For each updated model, the following "z statistic" is computed
  \deqn{z(\theta) = (\theta - \hat{\theta}) *
    \sqrt{\frac{D_{theta} - D_{\hat{theta}}}{\delta}}}{
    z(theta) = (theta - theta.hat) * sqrt((D_theta -
    D_theta.hat)/delta)} 
  where \eqn{\theta}{theta} is the constrained value of the parameter;
  \eqn{\hat{\theta}}{theta.hat} is the original fitted value;
  \eqn{D_{\theta}}{D_theta} is the deviance when the parameter is equal
  to \eqn{\theta}{theta}, and \eqn{\delta}{delta} is the dispersion
  parameter.

  When the deviance is quadratic in \eqn{\theta}{theta}, z will be
  linear in \eqn{\theta}{theta}. Therefore departures from quadratic
  behaviour can easily be identified by plotting z against
  \eqn{\theta}{theta} using \code{plot.profile.gnm}.

  \code{confint.profile.gnm} estimates confidence intervals for the
  parameters by interpolating the deviance profiles and identifying the
  parameter values at which z is equal to the relevant percentage points
  of the normal distribution. The \code{alpha} argument to
  \code{profile.gnm} specifies the significance level of z which must be
  covered by the profile. In particular, the profiling in a given
  direction will stop when \code{maxsteps} is reached or two steps have
  been taken in which
  \deqn{z(\theta) > (\theta - \hat{\theta}) * z_{(1 - \alpha)/2}}{
    z(theta) > (theta - theta.hat) * z_{(1 - alpha)/2}}

  By default, the stepsize is
  \deqn{z_{(1 - \alpha)/2} * s_{\hat{\theta}}}{
    z_{(1 - alpha)/2} * s_theta.hat}
  where \eqn{s_{\hat{\theta}}}{s_theta.hat} is the standard error of
  \eqn{\hat{\theta}}{theta.hat}. Strong assymetry is detected and
  the stepsize is adjusted accordingly, to try to ensure that the range
  determined by \code{alpha} is adequately covered. \code{profile.gnm}
  will also attempt to detect if the deviance is asymptotic such that
  the desired significance level cannot be reached. Each profile has an
  attribute \code{"asymptote"}, a two-length logical vector specifying
  whether an asymptote has been detected in either direction.

  For unidentified parameters the profile will be \code{NA}, as such
  parameters cannot be profiled.  
}
\value{
  A list of profiles, with one named component for each parameter
  profiled. Each profile is a data.frame: the first column, "z", contains
  the z statistics and the second column "par.vals" contains a matrix of
  parameter values, with one column for each parameter in the model.

  The list has two attributes: "original.fit" containing \code{fitted}
  and "summary" containing \code{summary(fitted)}.
}
\references{
  Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S} }
\author{ Heather Turner }
\seealso{ \code{\link{confint.gnm}},  \code{\link{gnm}},
  \code{\link[MASS]{profile.glm}}, \code{\link{ofInterest}}
}
\examples{
set.seed(1)

### Example in which deviance is near quadratic
data(voting)
count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)
classMobility <- gnm(yvar ~ -1 + Dref(origin, destination),
                     constrain = "delta1", family = binomial,
                     data = voting)
prof <- profile(classMobility, trace = TRUE)
plot(prof)
## confint similar to MLE +/- 1.96*s.e. 
confint(prof, trace = TRUE)
coefData <- se(classMobility)
cbind(coefData[1] - 1.96 * coefData[2], coefData[1] + 1.96 * coefData[2])

\dontrun{
### These examples take longer to run
### Another near quadratic example
data(occupationalStatus)
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination),
               ofInterest = "MultHomog", constrain = "MultHomog.*1",
               family = poisson, data = occupationalStatus)
prof <- profile(RChomog, trace = TRUE)
plot(prof)
## confint similar to MLE +/- 1.96*s.e. 
confint(prof)
coefData <- se(RChomog)
cbind(coefData[1] - 1.96 * coefData[2], coefData[1] + 1.96 * coefData[2])

## Another near quadratic example, with more complex constraints
data(voting)
count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)
classMobility <- gnm(yvar ~ -1 + Dref(origin, destination), 
                     family = binomial, data = voting)
wts <- prop.table(exp(coef(classMobility))[1:2])
classMobility <- update(classMobility, constrain = "delta1",
                        constrainTo = log(wts[1]))
sum(exp(parameters(classMobility))[1:2]) #=1
prof <- profile(classMobility, trace = TRUE)
plot(prof)
## confint similar to MLE +/- 1.96*s.e. 
confint(prof, trace = TRUE)
coefData <- se(classMobility)
cbind(coefData[1] - 1.96 * coefData[2], coefData[1] + 1.96 * coefData[2])

### An example showing asymptotic deviance
data(yaish)
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(educ), orig:dest),
               ofInterest = "[.]educ", constrain = "[.]educ1",
               family = poisson, data = yaish, subset = (dest != 7))
prof <- profile(unidiff, trace = TRUE)
plot(prof)
## clearly not quadratic for Mult1.Factor1.educ4 or Mult1.Factor1.educ5!
confint(prof)
##                          2.5 \%     97.5 \%
## Mult(Exp(.), orig:dest).educ1         NA         NA
## Mult(Exp(.), orig:dest).educ2 -0.5978901  0.1022447
## Mult(Exp(.), orig:dest).educ3 -1.4836854 -0.2362378
## Mult(Exp(.), orig:dest).educ4 -2.5792398 -0.2953420
## Mult(Exp(.), orig:dest).educ5       -Inf -0.7006889
## coefData <- se(unidiff)
## cbind(coefData[1] - 1.96 * coefData[2], coefData[1] + 1.96 * coefData[2])

### A far from quadratic example, also with eliminated parameters
data(backPain)
backPainLong <- expandCategorical(backPain, "pain")

oneDimensional <- gnm(count ~ pain + Mult(pain, x1 + x2 + x3),
                      eliminate = id,  family = "poisson",
                      constrain = "[.](painworse|x1)", constrainTo = c(0, 1),
                      data = backPainLong)
prof <- profile(oneDimensional, trace = TRUE)
plot(prof)
## not quadratic for any non-eliminated parameter
confint(prof)
coefData <- se(oneDimensional)
cbind(coefData[1] - 1.96 * coefData[2], coefData[1] + 1.96 * coefData[2])
}
}
\keyword{ models }
\keyword{ nonlinear }
