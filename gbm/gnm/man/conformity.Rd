\name{conformity}
\alias{conformity}
\docType{data}
\title{Conformity to Parental Rules}
\description{
  Data from a study on the value parents place on their children
  conforming to their rules.
}
\usage{data(conformity)}
\format{
  A data frame with 589 observations on the following 10 variables.
  \describe{
    \item{MCFM}{mother's conformity score.}
    \item{FCFF}{father's conformity score.}
    \item{MOPLM}{a factor describing the mother's education with levels
      \code{1} (low) to \code{7} (high).}
    \item{FOPLF}{a factor describing the father's education with levels
      \code{1} (low) to \code{7} (high).}
    \item{AGEM}{mother's birth cohort.}
    \item{MRMM}{mother's traditional role model.}
    \item{FRMF}{father's traditional role model.}
    \item{MWORK}{mother's employment.}
    \item{MFCM}{mother's family conflict score.}
    \item{FFCF}{father's family conflict score.}
  }
}
\details{
  The covariates \code{AGEM}, \code{MRMM}, \code{FRMF}, \code{MWORK},
  \code{MFCM} and \code{FFCF} have been recoded: for each covariate,
  values below or equal to the mean were recoded as zero,
  values above the mean were recoded as one.
}
}
\references{
  Van der Slik, F. W. P., De Graaf, N. D and Gerris, J. R. M. (2002)
  Conformity to Parental Rules: Asymmetric Influences of Father's and
  Mother's Levels of Education. \emph{European Sociological Review}
  \bold{18(4)}, 489 -- 502. 
}
\examples{
### Examples from Van der Slik paper
data(conformity)

## create factors required for status maximisation model
status.max <- as.numeric(MOPLM) < as.numeric(FOPLF)
MOPLMadj <- ifelse(status.max, FOPLF, MOPLM)
FOPLFadj <- ifelse(status.max, MOPLM, FOPLF)

## Models for mothers' conformity score
A <- gnm(MCFM ~ Nonlin(Dref(MOPLM, FOPLF)) + MFCM + AGEM + FRMF +
         MRMM + MWORK - 1, family = gaussian, data = slik, trace = TRUE)
exp(coef(A)[6])/(exp(coef(A)[6]) + exp(coef(A)[7])) # as in Table 4

B <- gnm(MCFM ~ Nonlin(Dref(MOPLMadj, FOPLFadj)) + MFCM
         + AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)

C <- gnm(MCFM ~ Nonlin(Dref(MOPLM, FOPLF, formula = ~ 1 + AGEM)) + MFCM +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)

D <- gnm(MCFM ~ Nonlin(Dref(MOPLM, FOPLF, formula = ~ 1 + MRMM)) + MFCM +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)

E <- gnm(MCFM ~ Nonlin(Dref(MOPLM, FOPLF, formula = ~ 1 + MWORK)) + MFCM +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)

F <- gnm(MCFM ~ Nonlin(Dref(MOPLM, FOPLF, formula = ~ 1 + MFCM)) + MFCM +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)

\dontrun{
## equivalent models for fathers
A <- gnm(FCFF ~ Nonlin(Dref(FOPLF, MOPLM)) + FFCF + AGEM + FRMF +
         MRMM + MWORK - 1, family = gaussian, data = slik)
exp(coef(A)[6])/(exp(coef(A)[6]) + exp(coef(A)[7])) # as in Table 4

B <- gnm(FCFF ~ Nonlin(Dref(FOPLFadj, MOPLMadj)) + FFCF
         + AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)

C <- gnm(FCFF ~ Nonlin(Dref(FOPLF, MOPLM, formula = ~ 1 + AGEM)) + FFCF +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)

D <- gnm(FCFF ~ Nonlin(Dref(FOPLF, MOPLM, formula = ~ 1 + MRMM)) + FFCF +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)

E <- gnm(FCFF ~ Nonlin(Dref(FOPLF, MOPLM, formula = ~ 1 + MWORK)) + FFCF +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)

F <- gnm(FCFF ~ Nonlin(Dref(FOPLF, MOPLM, formula = ~ 1 + FFCF)) + FFCF +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = slik)
}
}
\keyword{datasets}
