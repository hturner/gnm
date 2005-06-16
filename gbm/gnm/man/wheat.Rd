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
mainEffects <- glm(yield.scaled ~ year + treatment, family = gaussian,
                   data = wheat)
bilinear1 <- gnm(yield.scaled ~ year + treatment + Mult(year, treatment),
                 family = gaussian, data = wheat)
anova(mainEffects, bilinear1)
\dontrun{
    ## The next two take a lot of iterations to converge...
    bilinear2 <- gnm(yield.scaled ~ year + treatment +
                     Mult(year, treatment, multiplicity = 2),
                     family = gaussian, data = wheat)
    bilinear3 <- gnm(yield.scaled ~ year + treatment +
                     Mult(year, treatment, multiplicity = 3),
                     family = gaussian, data = wheat)
    anova(mainEffects, bilinear1, bilinear2, bilinear3)
}

##  Examine the extent to which, say, mTF explains the first bilinear term
bilinear1mTF <- gnm(yield.scaled ~ year + treatment + Mult(1 + mTF, treatment),
                    family = gaussian, data = wheat)
anova(mainEffects, bilinear1mTF, bilinear1)
}
\keyword{datasets}
