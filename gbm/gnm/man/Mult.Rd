\name{Multiplicative interaction}
\alias{Mult}
\alias{Exp}
\title{Specify a Multiplicative Interaction in a Model Formula}
\description{
\code{Mult} is a symbolic wrapper to indicate a multiplicative interaction in
the formula argument to \code{gnm}.

\code{Exp} is a symbolic wrapper to indicate an exponeniated multiplier
in a \code{Mult} term. 
}
\usage{
Mult(..., multiplicity = 1)
Exp(multiplier)
}
\arguments{
  \item{\dots}{a comma-separated list of two or more symbolic
    expressions, which represent the constituent multipliers in each
    component of the interaction.}
  \item{multiplicity}{a positive integer specifying the number of
    components to fit.}
  \item{multiplier}{a symbolic expression.}
}
\details{
  A multiplicative interaction is a sum of components that are of the
  form 
  \deqn{m_1 m_2 ... m_n}
  where the constituent multipliers \eqn{m_1, m_2, ..., m_n} are linear
  predictors, which are of the same form in each component.

  \code{Mult} may be used to specify a multiplicative interaction in the
  \code{formula} argument to \code{gnm}. The desired number of
  components in the interaction is specified by \code{multiplicity}.

  Models for the constituent multipliers are specified symbolically
  as unspecified arguments to \code{Mult}. These symbolic expressions
  are interpreted in the same way as the right hand side of a formula in
  an object of class \code{"formula"}. In particular, an intercept term
  is added by default. Offsets can be added to constituent multipliers,
  using \code{offset}.

  \code{Exp} is a companion function to \code{Mult}, which may be used to
  specify an exponentiated constituent multiplier, so that the
  multiplier is constrained to be nonnegative. 
}
\value{
    \code{Mult} returns a list of objects of class \code{"Mult"}, with one element
    for each multiplicative component. A \code{"Mult"} object is a list of deparsed
    expressions representing the factors in the multiplicative interaction.

    \code{Exp} returns an object of class \code{"Exp"} which is a deparsed
    expression representing the factor to be exponentiated.
  }
\references{
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
## Using 'Mult' with 'Exp' to constrain the first constituent multiplier
data(yaish)
## Fit the "UNIDIFF" mobility model across education levels
unidiff <- gnm(Freq ~ educ:orig + educ:dest +
               Mult(Exp(-1 + educ), orig:dest),
               family = poisson, data = yaish)

\dontrun{
## (this example can take quite a while to run)
##
## Using 'multiplicity' > 1
data(wheat)
yield.scaled <- wheat$yield * sqrt(3/1000)
treatment <- factor(paste(wheat$tillage, wheat$summerCrop, wheat$manure,
                          wheat$N, sep = ""))
bilinear2 <- gnm(yield.scaled ~ year + treatment +
                 Mult(year, treatment, multiplicity = 2),
                 family = gaussian, data = wheat)
}
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
