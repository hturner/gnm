\name{term.predictors}
\alias{term.predictors}
\title{ Extract Term Contributions to Predictor }
\description{
  \code{term.predictors} is a generic function which extracts the
  contribution of each term to the predictor from a fitted model object.  
}
\usage{
term.predictors(object, ...)
}
\arguments{
  \item{object}{ a fitted model object. }
  \item{\dots}{ additional arguments for method functions.  }
}
\details{
  The default method assumes that the predictor is linear and calculates
  the contribution of each term from the model matrix and fitted
  coefficients. A method is also available for \code{gnm} objects.
}
\value{
  A matrix with the additive components of the predictor in labelled columns.
}
\author{ Heather Turner }

\seealso{ \code{\link{gnm}} }
\examples{
## Linear model
G <- gl(4, 6)
x <- 1:24
y <- rnorm(24, 0, 1)
lmGx <- lm(y ~ G + x)
contrib <- term.predictors(lmGx)
contrib
all.equal(as.numeric(rowSums(contrib)), as.numeric(lmGx$fitted)) #TRUE

## Generalized linear model
y <- cbind(rbinom(24, 10, 0.5), rep(10, 24))
glmGx <- glm(y ~ G + x, family = binomial)
contrib <- term.predictors(glmGx)
contrib
all.equal(as.numeric(rowSums(contrib)),
as.numeric(glmGx$linear.predictors)) #TRUE

## Generalized nonlinear model
A <- gl(4, 6)
B <- gl(6, 1, 24)
y <- cbind(rbinom(24, 10, 0.5), rep(10, 24))
\dontrun{
## Using 'gnm' package
library(gnm)
}
gnmAB <- gnm(y ~ A + B + Mult(A - 1, B - 1), family = binomial)
contrib <- term.predictors(gnmAB)
contrib
all.equal(as.numeric(rowSums(contrib)),
as.numeric(gnmAB$predictors)) #TRUE
}

\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
