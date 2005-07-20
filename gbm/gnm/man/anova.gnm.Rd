\name{anova.gnm}
\alias{anova.gnm}
\title{ Analysis of Deviance for Generalized Nonlinear Models }
\description{
  Compute an analysis of deviance table for one or more generalized
  nonlinear models
}
\usage{
anova.gnm(object, ..., dispersion = NULL, test = NULL)
}
\arguments{
  \item{object}{ an object of class \code{gnm} }
  \item{\dots}{ additional objects of class \code{gnm} }
  \item{dispersion}{ the dispersion parameter for the fitting family. By
  default it is derived from \code{object} }
  \item{test}{ (optional) a character string, (partially) matching one
    of "Chisq", "F", or "Cp". See \code{\link{stat.anova}}. }
}
\details{
  Specifying a single object gives a sequential analysis of deviance
  table for that fit. The rows of the table show the reduction
  in the residual deviance and the current residual deviance as each
  term in the formula is added in turn.

  If more than one object is specified, the rows of the table show the
  residual deviance of the current model and the change in the residual
  deviance from the previous model. (This only makes statistical
  sense if the models are nested.) It is conventional to list the models
  from smallest to largest, but this is up to the user.

  If \code{test} is specified, the table will include test statistics
  and/or p values for the reduction in deviance. For models with known
  dispersion (e.g., binomial and Poisson fits) the chi-squared test is
  most appropriate, and for those with dispersion estimated by moments
  (e.g., 'gaussian', 'quasibinomial' and 'quasipoisson' fits) the F test
  is most appropriate.  Mallows' Cp statistic is the residual deviance
  plus twice the estimate of sigma^2 times the residual degrees of
  freedom, which is closely related to AIC (and a multiple of it if the
  dispersion is known). 
}
\value{
  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}
\author{ Heather Turner }
\section{Warning }{ The comparison between two or more models will only
  be valid if they are fitted to the same dataset. This may be a problem
  if there are missing values and R's default of 'na.action = na.omit'
  is used; an error will be given in this case. } 

\seealso{ \code{\link{gnm}}, \code{\link{anova}}}
\examples{
set.seed(1)
data(occupationalStatus)

##  Fit a uniform association model separating diagonal effects
Rscore <- scale(as.numeric(row(occupationalStatus)), scale = FALSE)
Cscore <- scale(as.numeric(col(occupationalStatus)), scale = FALSE)
Uniform <- gnm(Freq ~ origin + destination + Diag(origin, destination) + 
               Rscore:Cscore, family = poisson, data = occupationalStatus)

anova(Uniform, test = "Chisq")

##  Fit an association model with homogeneous row-column effects
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               Nonlin(MultHomog(origin, destination)), family = poisson,
               data = occupationalStatus)

##  Fit an association model with separate row and column effects
RC <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
          Mult(origin, destination), family = poisson,
          data = occupationalStatus)

anova(Uniform, RChomog, RC, test = "Chisq")

}
\keyword{ models }
