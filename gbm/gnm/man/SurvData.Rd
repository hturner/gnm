\name{SurvData}
\alias{SurvData}
\title{ Fill in history for right-censored survival data }
\description{
  Given the event status and follow-up time, \code{SurvData} will fill
  in the history of these variables and any other covariates, including
  simple time-dependent variables, producing a data frame with one row
  for each unit of time. 
}
\usage{
SurvData(time, event, data, timeDependent = clockControl(...), ...)
or
SurvData(event, data =, timeDependent = clockControl(...), ...)
 }
\arguments{
  \item{time}{ the follow-up times, specified either by a character
    vector which is the name of a column in \code{data}, or by a numeric
    vector }
  \item{event}{ if \code{time} is not specified, a \code{\link{Surv}}
    object giving both event status and follow-up times; otherwise the
    status indicator, specified in the same way as \code{time}. Status
    should be indicated by one of the following codings for
    {no event}/{event}: 0/1; FALSE/TRUE, or 1/2. If \code{event} is not
    specified, all subjects are assumed to have had an event at
    \code{time}. }
  \item{data}{ a data frame of variables recorded for each subject,
    possibly including the event time and status.  }
  \item{timeDependent}{ a list of parameters specifying the
    time-dependence of variables in \code{data}, where relevant. See
    \code{\link{clockControl}} for details.  }
  \item{\dots}{ arguments to be passed on to \code{clockControl} }
}
\details{
  Any time-dependent variables, including the follow-up times, are
  rounded up to the nearest integer. A record is then produced for each 
  unit of time from 1 up to \code{time}, for each subject.

  The amplified time and status data will appear in the returned data
  frame with the names "time" and "status", or the names of the vectors
  passed to \code{time} and \code{event} where applicable. If
  \code{data} contains variables with these names that are not identical
  to the variables specified by \code{time} and \code{status}, unique
  names for the amplified time and status data will be obtained using
  \code{make.unique}. 
  
  The function tries to distinguish between the use of 0/1 and 1/2
  coding for using \code{if (max(status) == 2)}. If 1/2 coding is used
  and all the subjects are censored, it will guess wrong.  Use 0/1
  coding in this case. 
}
\value{
  A data frame of the amplified time and status data, along with the
  corresponding covariate data.
}
\author{ Heather Turner }
\seealso{ \code{\link{clockControl}}, \code{\link{Surv}}}
\examples{
library(survival)
aml
SurvData("time", "status", data = aml)
}
\keyword{ manip }
\keyword{ survival }
