\name{meanResiduals}
\alias{meanResiduals}
\title{Average Residuals within Factor Levels}
\description{
  Computes the mean working residuals from a model fitted using
  Iterative Weighted Least Squares for each level of a factor or
  interaction of factors.
}
\usage{meanResiduals(object, by, standardized=TRUE, as.table=TRUE, ...)}
\arguments{
  \item{object}{model object for which \code{object$residuals} gives the
    working residuals and \code{object$weights} gives the working weights.}
  \item{by}{either a formula specifying a factor or interaction of
    factors (recommended), or a list of factors (the elements of which
    must correspond exactly to observations in the model frame). When a
    list of factors is specified, their interaction is used to specify
    the grouping factor.}
  \item{standardized}{logical: if \code{TRUE}, the mean residuals are
    standardized to be approximately standard normal.}
  \item{as.table}{logical: logical: if \code{TRUE} and \code{by}
    specifies an interaction of factors, the result is returned as a
    table cross-classified by these factors.}
  \item{...}{currently ignored}
}
\details{
  For level \eqn{i} of the grouping factor \eqn{A} the mean working
  residual is defined as
  \deqn{\frac{r_{ij} * w_{ij}}{\sum_{j = 1}^{n_i} w_{ij}}}{
    (r_ij * w_ij)/(sum_(j = 1)^(n_i) w_ij)}
  where \eqn{r_{ij}}{r_ij} is the \eqn{j}'th residual for level \eqn{i},
  \eqn{w_{ij}}{w_ij} is the corresponding working weight and \eqn{n_i}
  is the number of observations for level \eqn{i}. The denominator gives
  the weight corresponding to mean residual.

  For non-aggregated residuals, i.e. when the factor has one level per
  observation, the residuals are the same as Pearson residuals.
}
\author{Heather Turner}
\value{An object of class \code{"meanResiduals"}, for which \code{print}
  and \code{summary} methods are provided. A \code{"meanResiduals"}
  object is a list containing the following elements:
  \item{ call }{ the call used to create the model object from which the
    mean residuals are derived. }
  \item{ by }{ a label for the grouping factor. }
  \item{ residuals }{ the mean residuals. }
  \item{ df }{ the degrees of freedom associated with the mean
    residuals. }
  \item{ standardized }{ the \code{standardized} argument. }
  \item{ weights }{ the weights corresponding to the mean residuals. }
}
\examples{
## Fit a conditional independence model, leaving out
## the uninformative subtable for dest == 7:
CImodel <- gnm(Freq ~ educ*orig + educ*dest, family = poisson,
               data = yaish, subset = (dest != 7))

## compute mean residuals over origin and destination
meanRes <- meanResiduals(CImodel, ~ orig:dest)
meanRes
summary(meanRes)

\dontrun{
## requires vcdExtra package
## display mean residuals for origin and destination
library(vcdExtra)
mosaic(CImodel, ~orig+dest)
}

## non-aggregated residuals
res1 <- meanResiduals(CImodel, ~ educ:orig:dest)
res2 <- residuals(CImodel, type = "pearson")
all.equal(as.numeric(res1), as.numeric(res2))
}

\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
