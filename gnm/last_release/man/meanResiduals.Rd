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
  \item{by}{either a formula specifying one or more factors (recommended)
    or a list of factors (the elements of which must correspond exactly to
    observations in the model frame). When more than one factor is specified,
    their interaction is used to specify the grouping factor.}
  \item{standardized}{logical: if \code{TRUE}, the mean residuals are
    standardized to be approximately standard normal.}
  \item{as.table}{logical: if \code{TRUE} the result is returned as a table
    cross-classified by the factors passed to \code{by}.}
  \item{...}{currently ignored}
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
}
\examples{data(yaish)
## Fit a conditional independence model, leaving out
## the uninformative subtable for dest == 7:
CImodel <- gnm(Freq ~ educ*orig + educ*dest, family = poisson,
               data = yaish, subset = (dest != 7))

## compute mean residuals over origin and destination
meanRes <- meanResiduals(CImodel, ~ orig + dest)
meanRes
summary(meanRes)

## non-aggregated residuals
res1 <- meanResiduals(CImodel, ~ educ + orig + dest)
res2 <- residuals(CImodel, type = "pearson")
all.equal(res1$residuals[,,], res2[,,])
}

\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
