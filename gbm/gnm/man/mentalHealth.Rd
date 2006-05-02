\name{mentalHealth}
\alias{mentalHealth}
\docType{data}
\title{ Data on Mental Health and Socioeconomic Status}
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
mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)
RC1model <- gnm(count ~ SES + MHS +
                Mult(-1 + SES, -1 + MHS),
                family = poisson, data = mentalHealth)
## Row scores and column scores are both unnormalized in this
## parameterization of the model 

## The scores can be normalized as in Agresti's eqn (9.15):
rowProbs <- with(mentalHealth, tapply(count, SES, sum) / sum(count))
colProbs <- with(mentalHealth, tapply(count, MHS, sum) / sum(count))
rowScores <- coef(RC1model)[pickCoef(RC1model, "Mult1.*SES")]
colScores <- coef(RC1model)[pickCoef(RC1model, "Mult1.*MHS")]
rowScores <- rowScores - sum(rowScores * rowProbs)
colScores <- colScores - sum(colScores * colProbs)
beta1 <- sqrt(sum(rowScores^2 * rowProbs))
beta2 <- sqrt(sum(colScores^2 * colProbs))
assoc <- list(beta = beta1 * beta2,
              mu = rowScores / beta1,
              nu = colScores / beta2)
}
\keyword{datasets}
