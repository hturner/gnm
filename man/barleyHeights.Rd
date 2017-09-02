\name{barleyHeights}
\alias{barleyHeights}
\docType{data}
\title{ Heights of Barley Plants }
\description{
  Average heights for 15 genotypes of barley recorded over 9 years.
}
\usage{barleyHeights}
\format{
  A data frame with 135 observations on the following 3 variables.
  \describe{
    \item{\code{height}}{average height over 4 replicates (cm)}
    \item{\code{year}}{a factor with 9 levels \code{1974} to \code{1982}}
    \item{\code{genotype}}{a factor with 15 levels \code{1:15}}
  }
}
\source{
  Aastveit, A. H. \& Martens, H. (1986). ANOVA interactions interpreted
  by partial least squares regression. \emph{Biometrics}, \bold{42}, 829--844.
}
\references{
  Chadoeuf, J \& Denis, J B (1991). Asymptotic variances for the
  multiplicative interaction model. \emph{J. App. Stat.} \bold{18(3)}, 331--353.
}
\examples{
set.seed(1)
## Fit AMMI-1 model
barleyModel <- gnm(height ~ year + genotype + Mult(year, genotype),
                   data = barleyHeights)

## Get row and column scores with se's
gamma <- getContrasts(barleyModel, pickCoef(barleyModel, "[.]y"),
                      ref = "mean", scaleWeights = "unit")
delta <- getContrasts(barleyModel, pickCoef(barleyModel, "[.]g"),
                      ref = "mean", scaleWeights = "unit")

## Corresponding CI's similar to Chadoeuf & Denis (1991) Table 8
## (allowing for change in sign)
gamma[[2]][,1] + (gamma[[2]][,2]) \%o\% c(-1.96, 1.96)
delta[[2]][,1] + (delta[[2]][,2]) \%o\% c(-1.96, 1.96)

## Multiplier of row and column scores
height <- matrix(scale(barleyHeights$height, scale = FALSE), 15, 9)
R <- height - outer(rowMeans(height), colMeans(height), "+")
svd(R)$d[1]
}
\keyword{datasets}
