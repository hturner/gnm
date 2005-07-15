\name{SurvData}
\alias{SurvData}
\title{ Fill in history for right-censored survival data }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SurvData(time, event, data, timeDependent = clockControl(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{ ~~Describe \code{time} here~~ }
  \item{event}{ ~~Describe \code{event} here~~ }
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{timeDependent}{ ~~Describe \code{timeDependent} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(time, event, data,
                     timeDependent = clockControl(...), ...) {
    argList <- as.list(match.call())
    
    if (missing(event)) {
        if (missing(time))
            stop("No event data!")
        else if (inherits(time, "Surv"))
            event <- time
        else
            event <- 1
    }

    colNames <- c("time", "status")
   
    if (inherits(event, "Surv")) {
        if (ncol(event) != 2)
            stop("'SurvData' can not be used on interval event data")
        time <- event$time
        event <- event$status
    }
    else {
        if (is.character(time)) {
            colNames[1] <- time
            time <- data[[time]]
        }
        else {
            if (!is.numeric(time))
                stop("'time' is not of the correct format")
            if (is.name(argList$time))
                colNames[1] <- deparse(argList$time)[1]
        }
        if (is.character(event)) {
            colNames[2] <- event
            event <- data[[event]]
        }
        else {
            if (is.logical(event))
                event <- as.numeric(event)
            else if (!is.numeric(event))
                stop("'event' is not of the correct format")
            if (max(event) == 2)
                event <- event - 1
            if(is.name(argList$event))
                colNames[2] <- deparse(argList$event)[1]
        }
    }
                    
    roundTime <- ceiling(time)
    timeHistory <- sequence(roundTime)
    eventHistory <- replace(numeric(length(timeHistory)), cumsum(roundTime),
                            event)

    if (missing(data)) {
        survData <- cbind(timeHistory, eventHistory)
        colnames(survData) <- colNames
    }
    else {
        survData <- data[rep(1:nrow(data), roundTime), ]
        if (!identical(data[[colNames[1]]], time)) {
            survData[rev(make.unique(c(colnames(data), colNames[1])))[1]] <-
                timeHistory
        }
        else
            survData[, colNames[1]] <- timeHistory
        if (!identical(data[[colNames[2]]], event)) {
            survData[rev(make.unique(c(colnames(data), colNames[2])))[1]] <-
                eventHistory
        }
        else
            survData[, colNames[2]] <- eventHistory
        if (length(timeDependent$ind)) { 
            survData[, timeDependent$ind] <-
                mapply(function(name, x, time, decreasing, by) {
                    x <- as.matrix(x)
                    if (decreasing)
                        unlist(mapply(seq, x + by * (time - 1), x,
                                      by = -abs(by)))
                    else
                        unlist(mapply(seq, x - by * (time - 1), x,
                                      by = abs(by)))
                },
                       name = colnames(data[,timeDependent$ind, drop = FALSE]),
                       x = ceiling(data[,timeDependent$ind, drop = FALSE]),
                       decreasing = timeDependent$decreasing,
                       by = timeDependent$by,
                       MoreArgs = list(time = roundTime))
        }
    }
    survData
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
