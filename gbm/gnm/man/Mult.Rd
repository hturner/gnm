\name{Multiplicative interaction}
\alias{Mult}
\alias{Exp}
\title{Specify a multiplicative interaction in a model formula}
\description{
\code{Mult} is a symbolic wrapper to indicate a multiplicative interaction in
the formula argument to \code{gnm}.

\code{Exp} is a symbolic wrapper to indicate an exponeniated factor in a
\code{Mult} term.
}
\usage{
Mult(..., multiplicity = 1)
Exp(factor)
}
\arguments{
  \item{\dots}{a comma-separated list of two or more symbolic
    expressions representing the factors in the interaction.}
  \item{multiplicity}{a positive integer specifying the number of
    multiplicative components to fit.}
  \item{factor}{a symbolic expression representing a factor in a
    multiplicative interaction.}
}
\details{
    The factors in a \code{Mult} term are represented by symbolic
    expressions that are interpreted in the same way as the right hand
    side of a formula in an object of class \code{"formula"}. In
    particular, an intercept term is added by default. Offsets can be added to
    individual factors, using \code{offset}. Factors can be constrained
    to be nonnegative through the use of \code{Exp}.
}
\value{
    \code{Mult} returns a list of objects of class \code{"Mult"}, with one element
    for each multiplicative component. A \code{"Mult"} object is a list of deparsed
    expressions representing the factors in the multiplicative interaction.

    \code{Exp} returns an object of class \code{"Exp"} which is a deparsed
    expression representing the factor to be exponentiated.
  }
\references{
  Yaish, M (1998). Opportunities, Little Change. Class Mobility in
  Israeli Society: 1974-1991.  D.Phil. Thesis, Nuffield College,
  University of Oxford.

  Vargas, M, Crossa, J, van Eeuwijk, F, Sayre, K D and Reynolds, M P
  (2001).  Interpreting treatment by environment interaction in agronomy
  trials.  \emph{Agronomy Journal} \bold{93}, 949--960.
}
\author{Heather Turner}
\note{ These functions should only be used in formulae passed to \code{gnm}. }
\seealso{\code{\link{offset}}, \code{\link{gnm}}, \code{\link{Dref}},
         \code{\link{MultHomog}}}
\examples{
set.seed(1)
## Using 'Mult' with 'Exp' to constrain first factor
data(yaish)
## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ:orig + educ:dest +
               Mult(Exp(-1 + educ), orig:dest),
               family = poisson, data = yaish)

## Using 'multiplicity' > 1
data(wheat)
yield.scaled <- wheat$yield * sqrt(3/1000)
treatment <- factor(paste(wheat$tillage, wheat$summerCrop, wheat$manure,
                          wheat$N, sep = ""))
bilinear2 <- gnm(yield.scaled ~ year + treatment +
                 Mult(year, treatment, multiplicity = 2),
                 family = gaussian, data = wheat)
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
