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

\seealso{\code{\link{gnm}}}
\examples{


}
\keyword{ models }
\keyword{ regression }
