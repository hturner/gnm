\name{mentalHealth}
\alias{mentalHealth}
\docType{data}
\title{ Mental health and socioeconomic status}
\description{
  A 2-way contingency table from a sample of residents of Manhattan.
  Classifying variables are child's mental impairment (\code{MHS}) and
  parents' socioeconomic status (\code{SES}).
}
\usage{data(mentalHealth)}
\format{
  A data frame with 24 observations on the following 3 variables.
  \describe{
    \item{\code{count}}{a numeric vector}
    \item{\code{SES}}{an ordered factor with levels \code{A} < \code{B} < \code{C} < \code{D} < \code{E} < \code{F}}
    \item{\code{MHS}}{an ordered factor with levels \code{well} < \code{mild} < \code{moderate} < \code{impaired}}
  }
}
\source{
From Agresti (2002, p381); originally in Srole et al. (1978, p289).
}
\references{
Agresti, A. (2002).  \emph{Categorical Data Analysis} (2nd edn).  New
York: Wiley.

Srole, L, Langner, T. S., Michael, S. T., Opler, M. K. and Rennie,
T. A. C. (1978), \emph{Mental Health in the Metropolis: The Midtown
  Manhattan Study}.  New York: NYU Press.
}
\examples{
set.seed(1)
data(mentalHealth)

##  Goodman Row-Column association model fits well (deviance 3.57, df 8)
mh.SES <- factor(mentalHealth$SES, ordered = FALSE)
mh.MHS <- factor(mentalHealth$MHS, ordered = FALSE)
mh.y <- mentalHealth$count
mh.RC1model <- gnm(mh.y ~ mh.SES + mh.MHS +
                   Mult(-1 + mh.SES, -1 + mh.MHS),
                   family = poisson)
## Row scores are parameters coef(mh.RC1model)[10:15]
## Column scores are coef(mh.RC1model)[16:19]
## -- both unnormalized in this parameterization of the model

## The scores can be normalized as in Agresti's eqn (9.15):
mh.rowProbs <- tapply(mh.y, mh.SES, sum) / sum(mh.y)
mh.colProbs <- tapply(mh.y, mh.MHS, sum) / sum(mh.y)
mh.rowScores <- coef(mh.RC1model)[10:15]
mh.colScores <- coef(mh.RC1model)[16:19]
mh.rowScores <- mh.rowScores - sum(mh.rowScores * mh.rowProbs)
mh.colScores <- mh.colScores - sum(mh.colScores * mh.colProbs)
mh.beta1 <- sqrt(sum(mh.rowScores^2 * mh.rowProbs))
mh.beta2 <- sqrt(sum(mh.colScores^2 * mh.colProbs))
mh.assoc <- list(beta = mh.beta1 * mh.beta2,
                 mu = mh.rowScores / mh.beta1,
                 nu = mh.colScores / mh.beta2)
}
\keyword{datasets}
