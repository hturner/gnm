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
    \item{\code{MHS}}{an ordered factor with levels \code{Well} < \code{Mild} < \code{Moderate} < \code{Impaired}}
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
data(mentalHealth)
##  Goodman Row-Column association model fits well (deviance 3.57, df 8)
SES <- factor(mentalHealth$SES, ordered = FALSE)
MHS <- factor(mentalHealth$MHS, ordered = FALSE)
y <- mentalHealth$count
RC1model <- gnm(y ~ SES + MHS + Mult(-1 + SES, -1 + MHS),
    family = poisson)
## Row scores are parameters coef(RC1model)[10:15]
## Column scores are coef(RC1model)[16:19]
## -- both unnormalized in this parameterization of the model
}
\keyword{datasets}
