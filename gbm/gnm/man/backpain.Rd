\name{backpain}
\alias{backpain}
\docType{data}
\title{ Back Pain Prognosis }
\description{
  Data from a study of patients suffering from back pain. Prognostic
  variables were recorded at presentation and progress was categorised
  three weeks after treatment.
}
\usage{data(backpain)}
\format{
  A data frame with 101 observations on the following 4 variables.
  \describe{
    \item{x1}{length of previous attack.}
    \item{x2}{pain change.}
    \item{x3}{lordosis.}
    \item{pain}{an ordered factor describing the progress of each
      patient with levels \code{worse} < \code{same} <
      \code{slight.improvement} < \code{moderate.improvement} <
      \code{marked.improvement} < \code{complete.relief}. }
  }
}
\source{
  http://ideas.repec.org/c/boc/bocode/s419001.html
}
\references{
  Anderson, J. A. (1984) Regression and Ordered Categorical
  Variables. \emph{J. R. Statist. Soc. B}, \bold{46(1)}, 1-30.
}
\examples{
set.seed(1)
data(backpain)

### Fit logistic family models described in Table 5 of Anderson (1984)

## Model progress by category only
no.relationship <- gnm(pain ~ pain, family = "multinomial", data =
                       backpain)
logLik(no.relationship)

## Model progress assuming common dependence on prognostic variables
one.dimensional <- gnm(pain ~ pain + Mult(pain - 1, x1 + x2 + x3 - 1), family =
                       "multinomial", data = backpain)
logLik(one.dimensional)

## Model progress assuming group-specific dependence on prognostic variables
three.dimensional <- gnm(pain ~ pain + pain:(x1 + x2 + x3), family =
                         "multinomial", data = backpain)
logLik(three.dimensional)
}
\keyword{datasets}
