Diag <- function(..., binary = FALSE){
    dots <- list(...)
    dots <- lapply(dots, as.factor)
    Levels <- sort(unique(unlist(lapply(dots, levels))))
    facMatrix <- sapply(dots, as.character)
    f <- function(row){
        if (all(is.na(row))) return(NA)
        if (all(!is.na(row)) && all(row == row[1])) return(row[1])
        row <- na.omit(row)
        if (!all(row == row[1])) return(".")
        return(NA)
    }
    if (inherits(facMatrix, "matrix"))
        result <- factor(apply(facMatrix, 1, f), levels = c(".", Levels))
    else
        result <- factor(f(facMatrix))
    if (binary) result <- ifelse(result == ".", 0, 1)
    result
}
