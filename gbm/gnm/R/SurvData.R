SurvData <- function(time, status, data, timeDependent = clockControl(...),
                     ...) {
    roundTime <- ceiling(data$time)
    timeHistory <- sequence(roundTime)
    statusHistory <- replace(numeric(length(timeHistory)),
                             cumsum(roundTime), data$status)
    if (length(timeDependent$ind)) {
        caseData <- data[rep(1:nrow(data), roundTime), -timeDependent$ind]
        caseData[, c(time, status)] <- data.frame(timeHistory, statusHistory)
        historicalData <-
            mapply(function(x, time, decreasing, by) {
                x <- as.matrix(x)
                if (decreasing)
                    unlist(mapply(seq, x + by * (time - 1), x, by = -abs(by)))
                else
                    unlist(mapply(seq, x - by * (time - 1), x, by = abs(by)))
            },
                   x = ceiling(data[,timeDependent$ind, drop = FALSE]),
                   decreasing = timeDependent$decreasing,
                   by = timeDependent$by,
                   MoreArgs = list(time = roundTime))
        data <- data.frame(caseData, historicalData)
    }
    else {
        data <- data[rep(1:nrow(data), roundTime),]
        data[, c(time, status)] <- data.frame(timeHistory, statusHistory)
    }
    data
}
    
    
