\name{wedderburn}
\alias{wedderburn}
\title{ Wedderburn quasi-variance family }
\description{
  Creates a \code{link{family}} object for use with \code{\link{glm}},
  \code{\link{gnm}}, etc., for the variance function 
  \eqn{[\mu(1-\mu)]^2} introduced by Wedderburn (1974) for response values in
  [0,1].
}
\usage{
wedderburn(link = "logit")
}
\arguments{
  \item{link}{ The name of a link function.  Allowed are "logit", "probit"
  and "cloglog". }
}
\value{
  An object of class \code{\link{family}}.
}
\references{ 
Gabriel, K R (1998).  Generalised bilinear regression.  \emph{Biometrika} 
\bold{85}, 689--700.

McCullagh, P and Nelder, J A (1989).  \emph{Generalized Linear Models} (2nd ed).  
Chapman and Hall.

Wedderburn, R W M (1974).  Quasilikelihood functions, generalized linear models 
and the Gauss-Newton method.  \emph{Biometrika} \bold{61}, 439--47.
  
 }
\author{ David Firth }
\note{   The reported deviance involves an arbitrary constant (see McCullagh and
  Nelder, 1989, p330); for estimating dispersion, use the Pearson chi-squared
  statistic instead. }

\seealso{ \code{\link{glm}}, \code{\link{family}} }
\examples{
data(barley)  ##  data from Wedderburn (1974), see ?barley
##  Fit Wedderburn's logit model with variance proportional to the
##  square of mu(1-mu)
barley.model <- glm(y ~ site + variety, family = wedderburn, data = barley)
fit <- fitted(barley.model)
print(sum((barley$y - fit)^2 / (fit * (1-fit))^2))
##  Agrees with the chi-squared value reported in McCullagh and Nelder 
##  (1989, p331), which differs slightly from Wedderburn's reported value.
##
##  Fit the biplot model as in Gabriel (1998, p694)
barley.rank2 <- gnm(y ~ -1 + Mult(site, variety, multiplicity = 2),
    family = wedderburn, data = barley)
barley.svd <- svd(matrix(barley.rank2$predictors, 10, 9))
barley.A <- sweep(barley.svd$v, 2, sqrt(barley.svd$d), "*")[, 1:2]
barley.B <- sweep(barley.svd$u, 2, sqrt(barley.svd$d), "*")[, 1:2]
##  These are essentially A and B as in Gabriel (1998, p694), from which
##  the biplot is made by
plot(rbind(barley.A, barley.B), pch = c(LETTERS[1:9], as.character(1:9), "X"))
##
##  Fit the double-additive model as in Gabriel (1998, p697)
variety.binary <- factor(match(barley$variety, c(2,3,6), nomatch = 0) > 0,
    labels = c("Rest", "2,3,6"))
barley.doubleAdditive <- gnm(y ~ variety + Mult(site, variety.binary),
    family = wedderburn, data = barley)
}
\keyword{ models }
