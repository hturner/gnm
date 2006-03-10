Symm <- function(...){
    dots <- list(...)
    if (any(diff(sapply(dots, length)) != 0)) stop(
                "arguments to symm() must all have same length")
    dots <- lapply(dots, as.factor)
    facMatrix <- sapply(dots, as.character)
    result <- factor(apply(facMatrix, 1,
                           function(row){
                               string <- paste(sort(row), collapse = "")
                               if (any(is.na(row))) is.na(string) <- TRUE
                               string
                           }
                           )
                     )
    result
}
