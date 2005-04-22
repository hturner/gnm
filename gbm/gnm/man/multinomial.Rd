\name{multinomial}
\alias{multinomial}
\title{ A dummy multinomial family function for use with gnm models only }
\description{
  A function which returns  \code{list(family = "multinomial")}
}
\usage{
multinomial()
}
\details{
  This is used in calls to \code{gnm} where the response variable is
  categorical and a multinomial log-linear model is to be fitted.  The
  actual model fitting is done via the standard Poisson trick.
}
\value{
  A list, with one component:
  \item{family}{\code{"multinomial"}}
}
\author{ David Firth }
\references{
  Anderson, J. A. (1984) Regression and Ordered Categorical
  Variables. \emph{J. R. Statist. Soc. B}, \bold{46(1)}, 1-30.
}
\seealso{\code{\link{gnm}}}
\examples{
data(backpain)

## Fit the stereotype model described in Table 5 of Anderson (1984)
set.seed(1)
one.dimensional.A <- gnm(pain ~ pain + Mult(pain - 1, x1 + x2 + x3 - 1),
                         family = "multinomial", data = backpain,
                         model = TRUE, startit = 3)

## Fit the model in the same way as 'gnm': reformulated as poisson model
set.seed(1)
one.dimensional.B <- gnm(.counts ~ pain + Mult(pain - 1, x1 + x2 + x3 - 1),
                         eliminate = ~ .rowID, family = "poisson",
                         data = one.dimensional.A$model, startit = 3)
}
\keyword{ models }
\keyword{ regression }
