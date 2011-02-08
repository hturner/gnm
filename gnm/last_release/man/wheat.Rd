\name{wheat}
\alias{wheat}
\docType{data}
\title{ Wheat Yields from  Mexican Field Trials }
\description{
Data from a 10-year experiment at the CIMMYT experimental station
located in the Yaqui Valley near Ciudad Obregon, Sonora, Mexico --- factorial
design using 24 treatments in all.  In each of the 10 years the experiment was 
arranged in a randomized complete block design with three replicates.
}
\usage{data(wheat)}
\format{
  A data frame with 240 observations on the following 33 variables.
  \describe{
    \item{yield}{numeric, mean yield in kg/ha for 3 replicates}
    \item{year}{a factor with levels \code{1988:1997}}
    \item{tillage}{a factor with levels \code{T} \code{t}}
    \item{summerCrop}{a factor with levels \code{S} \code{s}}
    \item{manure}{a factor with levels \code{M} \code{m}}
    \item{N}{a factor with levels \code{0} \code{N} \code{n}}
    \item{MTD}{numeric, mean max temp sheltered (deg C) in December}
    \item{MTJ}{same for January}
    \item{MTF}{same for February}
    \item{MTM}{same for March}
    \item{MTA}{same for April}
    \item{mTD}{numeric, mean min temp sheltered (deg C) in December}
    \item{mTJ}{same for January}
    \item{mTF}{same for February}
    \item{mTM}{same for March}
    \item{mTA}{same for April}
    \item{mTUD}{numeric, mean min temp unsheltered (deg C)in December}
    \item{mTUJ}{same for January}
    \item{mTUF}{same for February}
    \item{mTUM}{same for March}
    \item{mTUA}{same for April}
    \item{PRD}{numeric, total precipitation (mm) in December}
    \item{PRJ}{same for January}
    \item{PRF}{same for February}
    \item{PRM}{same for March}
    \item{SHD}{numeric, mean sun hours in December}
    \item{SHJ}{same for January}
    \item{SHF}{same for February}
    \item{EVD}{numeric, total evaporation (mm) in December}
    \item{EVJ}{same for January}
    \item{EVF}{same for February}
    \item{EVM}{same for March}
    \item{EVA}{same for April}
  }
}
\source{ Tables A1 and A3 of
Vargas, M, Crossa, J, van Eeuwijk, F, Sayre, K D and Reynolds, M P
(2001).  Interpreting treatment by environment interaction in agronomy
trials.  \emph{Agronomy Journal} \bold{93}, 949--960.
}
\examples{
set.seed(1)
data(wheat)

##  Scale yields to reproduce analyses reported in Vargas et al (2001)
yield.scaled <- wheat$yield * sqrt(3/1000)

##  Reproduce (up to error caused by rounding) Table 1 of Vargas et al (2001)
aov(yield.scaled ~ year*tillage*summerCrop*manure*N, data = wheat)
treatment <- interaction(wheat$tillage, wheat$summerCrop, wheat$manure,
                         wheat$N, sep = "")
mainEffects <- lm(yield.scaled ~ year + treatment, data = wheat)
svdStart <- residSVD(mainEffects, year, treatment, 3)
bilinear1 <- update(asGnm(mainEffects), . ~ . + 
                    Mult(year, treatment),
                    start = c(coef(mainEffects), svdStart[,1]))
bilinear2 <- update(bilinear1, . ~ . + 
                    Mult(year, treatment, inst = 2),
                    start = c(coef(bilinear1), svdStart[,2]))
bilinear3 <- update(bilinear2, . ~ . + 
                    Mult(year, treatment, inst = 3),
                    start = c(coef(bilinear2), svdStart[,3]))
anova(mainEffects, bilinear1, bilinear2, bilinear3)

##  Examine the extent to which, say, mTF explains the first bilinear term
bilinear1mTF <- gnm(yield.scaled ~ year + treatment + Mult(1 + mTF, treatment),
                    family = gaussian, data = wheat)
anova(mainEffects, bilinear1mTF, bilinear1)

##  How to get the standard SVD representation of an AMMI-n model
##
##  We'll work with the AMMI-2 model, which here is called "bilinear2"
##
##  First, extract the contributions of the 5 terms in the model:
##
wheat.terms <- termPredictors(bilinear2)
##
##  That's a matrix, whose 4th and 5th columns are the interaction terms
##
##  Combine those two interaction terms, to get the total estimated
##  interaction effect:
##
wheat.interaction <- wheat.terms[, 4] + wheat.terms[, 5]
##
##  That's a vector, so we need to re-shape it as a 24 by 10 matrix
##  ready for calculating the SVD:
##
wheat.interaction <- matrix(wheat.interaction, 24, 10)
##
##  Now we can compute the SVD:
##
wheat.interaction.SVD <- svd(wheat.interaction)
##
##  Only the first two singular values are nonzero, as expected
##  (since this is an AMMI-2 model, the interaction has rank 2)
##
##  So the result object can be simplified by re-calculating the SVD with
##  just two dimensions:
##
wheat.interaction.SVD <- svd(wheat.interaction, nu = 2, nv = 2)
}
\keyword{datasets}
