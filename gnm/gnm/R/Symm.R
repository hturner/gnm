Symm <- function(...){
    dots <- list(...)
    if (any(diff(sapply(dots, length)) != 0)) stop(
                "arguments to symm() must all have same length")
    dots <- lapply(dots, as.factor)
    facMatrix <- sapply(dots, as.character)
    f <- function(row){
        string <- paste(sort(row), collapse = "")
        if (any(is.na(row))) is.na(string) <- TRUE
        string
    }
    if (inherits(facMatrix, "matrix"))
        result <- factor(apply(facMatrix, 1, f))
    else
        result <- factor(f(facMatrix))
    result
}
