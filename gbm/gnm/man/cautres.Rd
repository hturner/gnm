\name{cautres}
\alias{cautres}
\docType{data}
\title{ Class, religion and vote in France}
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
data(cautres)
##  Fit a "double UNIDIFF" model with the religion-vote and class-vote
##  interactions both modulated by nonnegative election-specific multipliers
doubleUnidiff <- gnm(Freq ~ election:vote + election:class:religion
    + Mult(Exp(election), religion:vote) + Mult(Exp(election), class:vote),
    family = poisson, data = cautres)
##
##  Examine the multipliers of the class-vote log odds ratios
##
coef(doubleUnidiff)[118:121]
#  Mult2.Factor1.election1 Mult2.Factor1.election2 Mult2.Factor1.election3 
#               -1.1899542              -0.5082191              -0.7405835 
#  Mult2.Factor1.election4 
#               -0.8281143
getContrasts(doubleUnidiff, 118:121)
#                            estimate        se
#  Mult2.Factor1.election1 -0.3618399 0.2534762
#  Mult2.Factor1.election2  0.3198951 0.1320034
#  Mult2.Factor1.election3  0.0875308 0.1446842
#  Mult2.Factor1.election4  0.0000000 0.0000000
##  Same thing but with election 1 as reference category:
getContrasts(doubleUnidiff, 121:118)
#                            estimate        se
#  Mult2.Factor1.election4  0.3618399 0.2534746
#  Mult2.Factor1.election3  0.4493707 0.2473524
#  Mult2.Factor1.election2  0.6817351 0.2401645
#  Mult2.Factor1.election1  0.0000000 0.0000000
}
\keyword{datasets}
