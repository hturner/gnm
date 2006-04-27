\name{cautres}
\alias{cautres}
\docType{data}
\title{ Data on Class, Religion and Vote in France}
\description{
  A 4-way contingency table of vote by class by religion
  in four French elections
}
\usage{data(cautres)}
\format{
    A table of counts, with classifying factors \code{vote} (levels
    \code{1:2}), \code{class} (levels \code{1:6}) and \code{religion}
    (levels \code{1:4}) and \code{election} (levels \code{1:4}).
}
\source{
  Bruno Cautres
}
\references{
Cautres, B, Heath, A F and Firth, D (1998).  Class, religion and vote in
Britain and France.  \emph{La Lettre de la Maison Francaise} \bold{8}.
}
\examples{
set.seed(1)
data(cautres)

##  Fit a "double UNIDIFF" model with the religion-vote and class-vote
##  interactions both modulated by nonnegative election-specific multipliers
doubleUnidiff <- gnm(Freq ~ election*vote + election*class*religion +
                     Mult(Exp(-1 + election), religion:vote) +
                     Mult(Exp(-1 + election), class:vote),
                     family = poisson, data = cautres)
##  Deviance should be 133.04

##  Examine the multipliers of the class-vote log odds ratios
coefs.of.interest <- grep("class.*\\.1", names(coef(doubleUnidiff)))
coef(doubleUnidiff)[coefs.of.interest]
##  Mult2.Factor1.election1 Mult2.Factor1.election2 
##               -0.5724370               0.1092972 
##  Mult2.Factor1.election3 Mult2.Factor1.election4 
##               -0.1230682              -0.2105843
##  Re-parameterize by setting Mult2.Factor1.election1 to zero
getContrasts(doubleUnidiff, coefs.of.interest)
## [[1]]
##                          estimate        SE    quasiSE    quasiVar
## Mult2.Factor1.election1 0.0000000 0.0000000 0.22854380 0.052232270
## Mult2.Factor1.election2 0.6817370 0.2401642 0.07395880 0.005469905
## Mult2.Factor1.election3 0.4493740 0.2473519 0.09475932 0.008979329
## Mult2.Factor1.election4 0.3618262 0.2534753 0.10934823 0.011957035

##  Same thing but with election 4 as reference category:
getContrasts(doubleUnidiff, rev(coefs.of.interest))
## [[1]]
##                            estimate        SE    quasiSE    quasiVar
## Mult2.Factor1.election4  0.00000000 0.0000000 0.10934823 0.011957035
## Mult2.Factor1.election3  0.08754785 0.1446834 0.09475932 0.008979329
## Mult2.Factor1.election2  0.31991082 0.1320023 0.07395880 0.005469905
## Mult2.Factor1.election1 -0.36182617 0.2534753 0.22854380 0.052232270
}
\keyword{datasets}
