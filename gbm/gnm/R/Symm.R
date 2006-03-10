Symm <- function(...){
    dots <- list(...)
    if (any(diff(sapply(dots, length)) != 0)) stop(
                "arguments to symm() must all have same length")
    dots <- lapply(dots, as.factor)
    Levels <- sort(unique(unlist(lapply(dots, levels))))
    facMatrix <- sapply(dots, as.character)
    nLevels <- length(Levels)
    combs <- gtools:::combinations(nLevels, length(dots), Levels,
                                   repeats.allowed = TRUE)
#    resultLevels <- apply(combs, 1, function(row) paste(row, collapse = ""))
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
