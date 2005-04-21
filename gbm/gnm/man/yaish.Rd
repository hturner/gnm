\name{yaish}
\alias{yaish}
\docType{data}
\title{ Class mobility by education in Israel}
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
\references{
  Yaish, M (1998). Opportunities, Little Change. Class Mobility in
  Israeli Society: 1974-1991.  D.Phil. Thesis, Nuffield College,
  University of Oxford. 
}
\examples{
set.seed(1)
data(yaish)

## Fit the "UNIDIFF" mobility model across education levels

yaish.unidiff <- gnm(Freq ~ educ:orig + educ:dest +
                     Mult(Exp(-1 + educ), orig:dest), family = poisson,
                     data = yaish)
## Deviance should be 208.3

##  Look at the multipliers of the orig:dest association:
coefs.of.interest <- grep("Mult1.*educ", names(coef(yaish.unidiff)))
coef(yaish.unidiff)[coefs.of.interest]
## Mult1.Factor1.educ1 Mult1.Factor1.educ2 Mult1.Factor1.educ3 
##           1.9853067           1.7685557           1.2540751 
## Mult1.Factor1.educ4 Mult1.Factor1.educ5 
##           0.9551895          -0.2525413

## Get standard errors for contrasts with educ5
getContrasts(yaish.unidiff, coefs.of.interest)
##                     estimate        se
## Mult1.Factor1.educ1 2.237848 0.9411054
## Mult1.Factor1.educ2 2.021097 0.9434948
## Mult1.Factor1.educ3 1.506616 0.9535574
## Mult1.Factor1.educ4 1.207731 0.9780784
## Mult1.Factor1.educ5 0.000000 0.0000000

## Table of model residuals
residuals(yaish.unidiff)
}
\author{David Firth}
\keyword{ datasets }

