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

## Fit the "UNIDIFF" mobility model across education levels

unidiff <- gnm(Freq ~ educ:orig + educ:dest +
                     Mult(Exp(-1 + educ), orig:dest), family = poisson,
                     data = yaish)
## Deviance should be 208.3

##  Look at the multipliers of the orig:dest association:
coefs.of.interest <- grep("Mult1.*educ", names(coef(unidiff)))
coef(unidiff)[coefs.of.interest]
## Mult1.Factor1.educ1 Mult1.Factor1.educ2 Mult1.Factor1.educ3 
##           1.9853067           1.7685557           1.2540751 
## Mult1.Factor1.educ4 Mult1.Factor1.educ5 
##           0.9551895          -0.2525413

## Get standard errors for contrasts with educ1
getContrasts(unidiff, coefs.of.interest)
## [[1]]
##                       estimate        SE   quasiSE    quasiVar
## Mult1.Factor1.educ1  0.0000000 0.0000000 0.0974357 0.009493714
## Mult1.Factor1.educ2 -0.2167488 0.1605291 0.1281367 0.016419027
## Mult1.Factor1.educ3 -0.7312150 0.2311479 0.2092979 0.043805616
## Mult1.Factor1.educ4 -1.0300838 0.3425260 0.3251657 0.105732710
## Mult1.Factor1.educ5 -2.2377782 0.9410555 0.9312462 0.867219422

## Table of model residuals
residuals(unidiff)
}
\author{David Firth}
\keyword{ datasets }

