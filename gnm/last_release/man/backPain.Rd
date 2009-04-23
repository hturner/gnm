\name{backPain}
\alias{backPain}
\docType{data}
\title{ Data on Back Pain Prognosis, from Anderson (1984) }
\description{
  Data from a study of patients suffering from back pain. Prognostic
  variables were recorded at presentation and progress was categorised
  three weeks after treatment.
}
\usage{data(backPain)}
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
  \url{http://ideas.repec.org/c/boc/bocode/s419001.html}
}
\references{
  Anderson, J. A. (1984) Regression and Ordered Categorical
  Variables. \emph{J. R. Statist. Soc. B}, \bold{46(1)}, 1-30.
}
\examples{
set.seed(1)
data(backPain)
summary(backPain)

### Re-express as count data
backPainLong <- expandCategorical(backPain, "pain")

### Fit models described in Table 5 of Anderson (1984)

### Logistic family models
noRelationship <- gnm(count ~ pain, eliminate = id,
                      family = "poisson", data = backPainLong)

## stereotype model
oneDimensional <- update(noRelationship,
                         ~ . + Mult(pain, x1 + x2 + x3))

## multinomial logistic
threeDimensional <- update(noRelationship, ~ . + pain:(x1 + x2 + x3))

### Models to determine distinguishability in stereotype model
## constrain scale of category-specific multipliers
oneDimensional <- update(noRelationship,
                         ~ . + Mult(pain, offset(x1) + x2 + x3))
## obtain identifiable contrasts & id possibly indistinguishable slopes
getContrasts(oneDimensional, pickCoef(oneDimensional, "[.]pain"))

\dontrun{
## (this part not needed for package testing)
## fit simpler models and compare 
.pain <- backPainLong$pain

levels(.pain)[2:3] <- paste(levels(.pain)[2:3], collapse = " | ")
fiveGroups <- update(noRelationship,
                     ~ . + Mult(.pain, x1 + x2 + x3))

levels(.pain)[4:5] <- paste(levels(.pain)[4:5], collapse = " | ")
fourGroups <- update(fiveGroups)

levels(.pain)[2:3] <- paste(levels(.pain)[2:3], collapse = " | ")
threeGroups <- update(fourGroups)

### Grouped continuous model, aka proportional odds model
library(MASS)
sixCategories <- polr(pain ~ x1 + x2 + x3, data = backPain)

### Obtain number of parameters and log-likelihoods for equivalent
### multinomial models as presented in Anderson (1984)
logLikMultinom <- function(model){
    object <- get(model)
    if (inherits(object, "gnm")) {
        l <- logLik(object) + object$eliminate
        c(nParameters = attr(l, "df") - object$eliminate, logLikelihood = l)
    }
    else
        c(nParameters = object$edf, logLikelihood = -deviance(object)/2)
}
models <- c("threeDimensional", "oneDimensional", "noRelationship",
            "fiveGroups", "fourGroups", "threeGroups", "sixCategories")
t(sapply(models, logLikMultinom))
}
}
\keyword{datasets}
