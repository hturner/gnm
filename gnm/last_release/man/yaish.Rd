\name{yaish}
\alias{yaish}
\docType{data}
\title{ Class Mobility by Level of Education in Israel}
\description{
  A 3-way contingency table of father/son pairs, classified by father's
  social class (\code{orig}), son's social class (\code{dest}) and son's 
  education level (\code{educ}).
}
\usage{data(yaish)}
\format{
A table of counts, with classifying factors \code{educ} (levels \code{1:5}),
\code{orig} (levels \code{1:7}) and \code{dest} (levels \code{1:7}).
}
\source{Originally in Yaish (1998), see also Yaish (2004, p316).}
\references{
  Yaish, M (1998). Opportunities, Little Change. Class Mobility in
  Israeli Society: 1974-1991.  D.Phil. Thesis, Nuffield College,
  University of Oxford. 
  
  Yaish, M (2004).  \emph{Class Mobility Trends in Israeli Society, 1974-1991.}
  Lewiston: Edwin Mellen Press.
}
\examples{
set.seed(1)

data(yaish)

## Fit the "UNIDIFF" mobility model across education levels, leaving out
## the uninformative subtable for dest == 7:
##
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
                     Mult(Exp(educ), orig:dest), family = poisson,
                     data = yaish, subset = (dest != 7))
## Deviance should be 200.3, 116 d.f.
##
## Look at the multipliers of the orig:dest association:
ofInterest(unidiff) <- pickCoef(unidiff, "[.]educ")
coef(unidiff)
##
## Coefficients of interest:
## Mult(Exp(.), orig:dest).educ1 Mult(Exp(.), orig:dest).educ2 
##                    -0.5513258                    -0.7766976 
## Mult(Exp(.), orig:dest).educ3 Mult(Exp(.), orig:dest).educ4 
##                    -1.2947494                    -1.5902644 
## Mult(Exp(.), orig:dest).educ5 
##                    -2.8008285 
##
## Get standard errors for the contrasts with educ1:
##
getContrasts(unidiff, ofInterest(unidiff))
##                                 estimate        SE    quasiSE
## Mult(Exp(.), orig:dest).educ1  0.0000000 0.0000000 0.09757438
## Mult(Exp(.), orig:dest).educ2 -0.2253718 0.1611874 0.12885847
## Mult(Exp(.), orig:dest).educ3 -0.7434236 0.2335083 0.21182123
## Mult(Exp(.), orig:dest).educ4 -1.0389386 0.3434256 0.32609380
## Mult(Exp(.), orig:dest).educ5 -2.2495026 0.9453764 0.93560643
##                                 quasiVar
## Mult(Exp(.), orig:dest).educ1 0.00952076
## Mult(Exp(.), orig:dest).educ2 0.01660450
## Mult(Exp(.), orig:dest).educ3 0.04486823
## Mult(Exp(.), orig:dest).educ4 0.10633716
## Mult(Exp(.), orig:dest).educ5 0.87535940
##
## Table of model residuals:
##
residuals(unidiff)
}
\author{David Firth}
\keyword{ datasets }

