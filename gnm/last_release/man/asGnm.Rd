\name{asGnm}
\alias{asGnm}
\title{ Coerce Linear Model to gnm Object }
\description{
  \code{asGnm} is a generic function which coerces objects of class
  "glm" or "lm" to an object of class "gnm". 
}
\usage{
asGnm(object, ...)
}
\arguments{
  \item{object}{ an object of class "glm" or "lm". }
  \item{\dots}{ additional arguments for method functions. }
}
\details{
  Components are added to or removed from \code{object} to produce an
  object of class "gnm". This can be useful in model building,
  see examples.
}
\value{
  An object of class "gnm" - see \code{\link{gnm}} for full description.
}
\references{
  Vargas, M, Crossa, J, van Eeuwijk, F, Sayre, K D and Reynolds, M P
  (2001).  Interpreting treatment by environment interaction in agronomy
  trials.  \emph{Agronomy Journal} \bold{93}, 949--960.
}
\author{ Heather Turner }
\seealso{ \code{\link{gnm}}, \code{\link{glm}}, \code{\link{lm}} }
\examples{
set.seed(1)

##  Scale yields to reproduce analyses reported in Vargas et al (2001)
yield.scaled <- wheat$yield * sqrt(3/1000)
treatment <- interaction(wheat$tillage, wheat$summerCrop, wheat$manure,
                         wheat$N, sep = "")
##  Fit linear model
mainEffects <- lm(yield.scaled ~ year + treatment, data = wheat)

##  Convert to gnm object to allow addition of Mult() term
svdStart <- residSVD(mainEffects, year, treatment, 3)
bilinear1 <- update(asGnm(mainEffects), . ~ . + 
                    Mult(year, treatment),
                    start = c(coef(mainEffects), svdStart[,1]))
}
\keyword{ models }
\keyword{ regression }
