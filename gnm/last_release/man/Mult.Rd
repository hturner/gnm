\name{Multiplicative interaction}
\alias{Mult}
\title{Specify a Product of Predictors in a gnm Model Formula}
\description{
A function of class \code{"nonlin"} to specify a multiplicative interaction in
the formula argument to \code{\link{gnm}}.

}
\usage{
Mult(..., inst = NULL)
}
\arguments{
  \item{\dots}{a comma-separated list of two or more symbolic
    expressions representing the constituent multipliers in the
    interaction.} 
  \item{inst}{a positive integer specifying the instance number of the
    term.} 
}
\details{
  \code{Mult} specifies instances of a multiplicative interaction,
  i.e. a product of the form 
  \deqn{m_1 m_2 ... m_n,}
  where the constituent multipliers \eqn{m_1, m_2, ..., m_n} are linear
  or nonlinear predictors.

  Models for the constituent multipliers are specified symbolically
  as unspecified arguments to \code{Mult}. These symbolic expressions
  are interpreted in the same way as the right hand side of a formula in
  an object of class \code{"formula"}, except that an intercept term
  is not added by default. Offsets can be added to constituent
  multipliers, using \code{offset}.

  The family of multiplicative interaction models include
  row-column association models for contingency tables (e.g., Agresti,
  2002, Sec 9.6), log-multiplicative or UNIDIFF models (Erikson and
  Goldthorpe, 1992; Xie, 1992), and GAMMI models (van Eeuwijk, 1995).
}
\value{
  A list with the required components of a \code{"nonlin"} function:
  \item{ predictors }{ the expressions passed to \code{Mult}}
  \item{ term }{ a function to create a deparsed mathematical expression
    of the term, given labels for the predictors.}
  \item{ call }{ the call to use as a prefix for parameter labels. }
}
\references{
  Agresti, A (2002).  \emph{Categorical Data Analysis} (2nd ed.)  New
  York: Wiley.
  
  Erikson, R and Goldthorpe, J H (1992).  \emph{The Constant Flux}.
  Oxford: Clarendon Press.

  van Eeuwijk, F A (1995).  Multiplicative interaction in generalized
  linear models.  \emph{Biometrics} \bold{51}, 1017-1032.

  Vargas, M, Crossa, J, van Eeuwijk, F, Sayre, K D and Reynolds, M P
  (2001).  Interpreting treatment by environment interaction in agronomy
  trials.  \emph{Agronomy Journal} \bold{93}, 949--960.

  Xie, Y (1992).  The log-multiplicative layer effect model for comparing
  mobility tables.  \emph{American Sociological Review} \bold{57},
  380-395.
}
\author{Heather Turner}
\seealso{\code{\link{gnm}}, \code{\link{formula}}, \code{\link{instances}},
  \code{\link{nonlin.function}}, \code{\link{MultHomog}}
}
\examples{
set.seed(1)

## Using 'Mult' with 'Exp' to constrain the first constituent multiplier
##  to be non-negative

## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(educ), orig:dest),
               family = poisson, data = yaish, subset = (dest != 7))

\dontrun{
## (this example can take quite a while to run)
##
## Fitting two instances of a multiplicative interaction (i.e. a
## two-component interaction))
yield.scaled <- wheat$yield * sqrt(3/1000)
treatment <- factor(paste(wheat$tillage, wheat$summerCrop, wheat$manure,
                          wheat$N, sep = ""))
bilinear2 <- gnm(yield.scaled ~ year + treatment +
                 instances(Mult(year, treatment), 2),
                 family = gaussian, data = wheat)
formula(bilinear2)
## yield.scaled ~ year + treatment + Mult(year, treatment, inst = 1) + 
##     Mult(year, treatment, inst = 2)
}
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
