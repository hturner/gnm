\name{barley}
\alias{barley}
\docType{data}
\title{ Jenkyn's Data on Leaf-blotch on Barley }
\description{
  Incidence of \emph{R. secalis} on the leaves of ten varieties of barley grown
  at nine sites.
}
\usage{data(barley)}
\format{
  A data frame with 90 observations on the following 3 variables.
  \describe{
    \item{y}{the proportion of leaf affected (values in [0,1])}
    \item{site}{a factor with 9 levels \code{A} to \code{I}}
    \item{variety}{a factor with 10 levels \code{c(1:9, "X")}}
  }
}
\note{
  This dataset was used in Wedderburn's original paper (1974) on 
  quasi-likelihood.
}
\source{
  Originally in an unpublished Aberystwyth PhD thesis by J F Jenkyn.
}
\references{
Gabriel, K R (1998).  Generalised bilinear regression.  \emph{Biometrika} 
\bold{85}, 689--700.

McCullagh, P and Nelder, J A (1989)  \emph{Generalized Linear Models} (2nd ed).  
Chapman and Hall.

Wedderburn, R W M (1974).  Quasilikelihood functions, generalized linear models 
and the Gauss-Newton method.  \emph{Biometrika} \bold{61}, 439--47.

}
\examples{
data(barley)
set.seed(1)

###  Fit Wedderburn's logit model with variance proportional to [mu(1-mu)]^2
logitModel <- glm(y ~ site + variety, family = wedderburn, data = barley)
fit <- fitted(logitModel)
print(sum((barley$y - fit)^2 / (fit * (1-fit))^2))
##  Agrees with the chi-squared value reported in McCullagh and Nelder 
##  (1989, p331), which differs slightly from Wedderburn's reported value.

###  Fit the biplot model as in Gabriel (1998, p694)
biplotModel <- gnm(y ~ -1 + Mult(site, variety, multiplicity = 2),
                   family = wedderburn, data = barley)
barleySVD <- svd(matrix(biplotModel$predictors, 10, 9))
A <- sweep(barleySVD$v, 2, sqrt(barleySVD$d), "*")[, 1:2]
B <- sweep(barleySVD$u, 2, sqrt(barleySVD$d), "*")[, 1:2]
##  These are essentially A and B as in Gabriel (1998, p694), from which
##  the biplot is made by
plot(rbind(A, B), pch = c(levels(barley$site), levels(barley$variety)))

##  Fit the double-additive model as in Gabriel (1998, p697)
variety.binary <- factor(match(barley$variety, c(2,3,6), nomatch = 0) > 0,
                        labels = c("rest", "2,3,6"))
doubleAdditive <- gnm(y ~ variety + Mult(site, variety.binary),
                      family = wedderburn, data = barley)
##  It is unclear why Gabriel's chi-squared statistics differ slightly
##  from the ones produced in these fits.  Possibly Gabriel adjusted the
##  data somehow prior to fitting?
}
\keyword{datasets}
