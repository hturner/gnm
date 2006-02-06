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
  
  Yaish, M (2004).  \em{Class Mobility Trends in Israeli Society, 1974-1991.}
  Lewiston: Edwin Mellen Press.
}
\examples{
set.seed(1)
data(yaish)

## Fit the "UNIDIFF" mobility model across education levels, leaving out
## the uninformative subtable for dest == 7:
##
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
                     Mult(Exp(-1 + educ), -1 + orig:dest), family = poisson,
                     data = yaish, subset = (dest != 7))
## Deviance should be 200.3, 116 d.f.
##
## Look at the multipliers of the orig:dest association:
coefs.of.interest <- grep("Mult1.*educ", names(coef(unidiff)))
coef(unidiff)[coefs.of.interest]
##
## Mult1.Factor1.educ1 Mult1.Factor1.educ2 Mult1.Factor1.educ3 Mult1.Factor1.educ4 
##          -0.8150195          -1.0403857          -1.5584007          -1.8538570 
## Mult1.Factor1.educ5 
##          -3.0643114
##
## Get standard errors for the contrasts with educ1:
##
getContrasts(unidiff, coefs.of.interest)
##
## [[1]]
##                       Estimate Std. Error   quasiSE   quasiVar
## Mult1.Factor1.educ1  0.0000000  0.0000000 0.0975751 0.00952090
## Mult1.Factor1.educ2 -0.2253661  0.1611886 0.1288595 0.01660478
## Mult1.Factor1.educ3 -0.7433812  0.2335042 0.2118164 0.04486617
## Mult1.Factor1.educ4 -1.0388374  0.3434049 0.3260719 0.10632288
## Mult1.Factor1.educ5 -2.2492918  0.9452228 0.9354512 0.87506893
##
## Table of model residuals:
##
residuals(unidiff)
}
\author{David Firth}
\keyword{ datasets }

