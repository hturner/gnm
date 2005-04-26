gnmControl <- function(tolerance = 1e-4, iterStart = 2, iterMax = 500,
                       trace = FALSE) {
    if (!is.numeric(iterMax) || iterMax <= 0) 
        stop("maximum number of iterations must be > 0")
    list(tolerance = tolerance, iterStart = iterStart, iterMax = iterMax,
         trace = trace)
}
