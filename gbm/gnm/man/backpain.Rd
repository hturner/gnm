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
data(backPain)

### Re-express as count data
.incidence <- class.ind(backPain$pain)
.counts <- as.vector(t(.incidence))
.rowID <- factor(t(row(.incidence)))
backPain <- backPain[.rowID, ]
backPain$pain <- factor(rep(levels(backPain$pain), nrow(.incidence)),
                        levels = levels(backPain$pain), ordered = TRUE)

### Fit models described in Table 5 of Anderson (1984)

### Logistic family models

noRelationship <- gnm(.counts ~ pain, eliminate = ~ .rowID,
                       family = "poisson", data = backPain)
logLik(noRelationship)

## stereotype model
oneDimensional <- update(noRelationship,
                         ~ . + Mult(pain - 1, x1 + x2 + x3 - 1),
                         startit = 3)
logLik(oneDimensional)

threeDimensional <- update(noRelationship, ~ . + pain:(x1 + x2 + x3))
logLik(threeDimensional)

### Models to determine distinguishability in stereotype model
pain2 <- backPain$pain

levels(pain2)[2:3] <- paste(levels(pain2)[2:3], collapse = " | ")
fiveGroups <- update(noRelationship,
                      ~ . + Mult(as.ordered(pain2) - 1,
                                 x1 + x2 + x3 - 1))
logLik(fiveGroups)

levels(pain2)[4:5] <- paste(levels(pain2)[4:5], collapse = " | ")
fourGroups <- update(fiveGroups)
logLik(fourGroups)

levels(pain2)[2:3] <- paste(levels(pain2)[2:3], collapse = " | ")
threeGroups <- update(fourGroups)
logLik(threeGroups)

### Grouped continuous model, aka proportional odds model
sixCategories <- polr(pain ~ x1 + x2 + x3, data = backpain)
-deviance(sixCategories)/2 # log-likelihood
}
\keyword{datasets}
