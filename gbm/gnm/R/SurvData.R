SurvData <- function(time, event, data,
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
    
    
