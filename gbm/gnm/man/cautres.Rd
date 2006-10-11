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
                     Mult(Exp(election), religion:vote) +
                     Mult(Exp(election), class:vote),
                     family = poisson, data = cautres)
##  Deviance should be 133.04

##  Examine the multipliers of the class-vote log odds ratios
ofInterest(doubleUnidiff) <- pickCoef(doubleUnidiff, "class:vote[).]")
coef(doubleUnidiff)
## Coefficients of interest:
## Mult(Exp(.), class:vote).election1 
##                        -0.38357138 
## Mult(Exp(.), class:vote).election2 
##                         0.29816599 
## Mult(Exp(.), class:vote).election3 
##                         0.06580307 
## Mult(Exp(.), class:vote).election4 
##                        -0.02174104
##  Re-parameterize by setting Mult2.Factor1.election1 to zero
getContrasts(doubleUnidiff, ofInterest(doubleUnidiff))
##                                     estimate        SE
## Mult(Exp(.), class:vote).election1 0.0000000 0.0000000
## Mult(Exp(.), class:vote).election2 0.6817374 0.2401644
## Mult(Exp(.), class:vote).election3 0.4493745 0.2473521
## Mult(Exp(.), class:vote).election4 0.3618301 0.2534754
##                                       quasiSE    quasiVar
## Mult(Exp(.), class:vote).election1 0.22854401 0.052232363
## Mult(Exp(.), class:vote).election2 0.07395886 0.005469913
## Mult(Exp(.), class:vote).election3 0.09475938 0.008979340
## Mult(Exp(.), class:vote).election4 0.10934798 0.011956981

##  Same thing but with election 4 as reference category:
getContrasts(doubleUnidiff, rev(ofInterest(doubleUnidiff)))
##                                       estimate        SE
## Mult(Exp(.), class:vote).election4  0.00000000 0.0000000
## Mult(Exp(.), class:vote).election3  0.08754436 0.1446833
## Mult(Exp(.), class:vote).election2  0.31990727 0.1320022
## Mult(Exp(.), class:vote).election1 -0.36183013 0.2534754
##                                       quasiSE    quasiVar
## Mult(Exp(.), class:vote).election4 0.10934798 0.011956981
## Mult(Exp(.), class:vote).election3 0.09475938 0.008979340
## Mult(Exp(.), class:vote).election2 0.07395886 0.005469913
## Mult(Exp(.), class:vote).election1 0.22854401 0.052232363
}
\keyword{datasets}
