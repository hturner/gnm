\name{clockControl}
\alias{clockControl}
\title{ Set Clock Parameters for Time-Dependent Variables }
\description{
  Set clock parameters for time-dependent variables in the data frame
  being passed to the \code{data} argument of \code{\link{SurvData}}.
}
\usage{
clockControl(ind = numeric(0), decreasing = FALSE, by = 1)
}
\arguments{
  \item{ind}{ either a numeric vector or character vector specifying the
    columns of the data frame which are time-dependent.  }
  \item{decreasing}{ a logical vector of length \code{1} or
    \code{length(ind)}, indicating whether the variables are decreasing
    over time.}
  \item{by}{ a numeric vector of length \code{1} or \code{length(ind)},
    specifying the increment of each variable from one time point to the
    next. }
}
\value{
  A list with the arguments as components
}
\author{ Heather Turner }
\seealso{ \code{\link{SurvData}} }
\examples{
Data <- data.frame(time = c(12.6, 3.1, 5.2, 10.1),
                   status = c(1, 0, 0, 1),
                   x = c("A", "F", "B", "E"),
                   age = c(28, 22, 35, 30))
SurvData("time", "status", Data, ind = 4)
}
\keyword{ manip }
\keyword{ survival }
