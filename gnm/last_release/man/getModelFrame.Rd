\name{getModelFrame}
\alias{getModelFrame}
\title{ Get the Model Frame in Use by 'gnm'   }
\description{
  A utility function to enable \code{gnm} plug-in functions to get the
  model frame in use by \code{\link{gnm}}.
}
\usage{
getModelFrame()
}
\value{
  The model frame in use by \code{gnm}.
}
\note{This function is only designed to work from within a plug-in
  function. To extract the model frame from a \code{gnm} object, use
  \code{\link{model.frame}}.}
\author{ Heather Turner }
\seealso{ \code{\link{gnm}}, \code{\link{Dref}}, \code{\link{MultHomog}} }
\examples{
## Create a dummy plug-in function
dummy <- function(...) {
    cat("Model frame:\n")
    print(getModelFrame()[1:5,])
    stop("Not a valid plug-in function, model can not be estimated")
}

##  Use data from Vargas et al (2001)
set.seed(1)
data(wheat)
yield.scaled <- wheat$yield * sqrt(3/1000)
treatment <- factor(paste(wheat$tillage, wheat$summerCrop, wheat$manure,
                          wheat$N, sep = ""))

##  Add dummy nonlinear term to main effects model - the dummy function
##  gets the model frame used by gnm and prints the first 5 rows
mainEffects <- try(gnm(yield.scaled ~ year + treatment +
                       Nonlin(dummy(as.numeric(year), N, scale(MTD))),
                       family = gaussian, data = wheat))
}
\keyword{ data }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
