gnm.control <- function(epsilon = 1e-4, startit = 10, maxit = 500,
                         trace = FALSE) {
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, startit = startit, maxit = maxit, trace = trace)
}
