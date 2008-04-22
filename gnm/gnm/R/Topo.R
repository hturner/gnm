"Topo" <-
    function (..., spec = NULL)
{
    if (is.null(spec)) stop("No spec given")
    dots <- list(...)
    factorLengths <- sapply(dots, length)
    lengthsEqual <- {if (length(factorLengths) == 1) TRUE else
                     sd(factorLengths) == 0}
    if (!lengthsEqual) stop("Factors have different lengths")
    specDim <- if (is.vector(spec)) length(spec) else dim(spec)
    dots <- lapply(dots, as.factor)
    facMat <- cbind(...)
    spec.ok <- identical(sapply(dots, nlevels), specDim)
    if (!spec.ok) stop(
            "Dimensions of spec do not match the factor arguments")
    return(as.factor(spec[facMat]))
}
