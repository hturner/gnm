#  Copyright (C) 2005, 2008 David Firth and Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

Diag <- function(..., binary = FALSE){
    dots <- list(...)
    dots <- lapply(dots, as.factor)
    Levels <- levels(dots[[1]])
    check <- vapply(dots[-1], function(x) identical(levels(x), Levels), 
                    TRUE)
    if (!all(check)){
        message("Levels are not identical, new factor will be based ",
                "on sorted combined levels.")
        Levels <- sort(unique(unlist(lapply(dots, levels))))
    }
    facMatrix <- vapply(dots, as.character, character(length(dots[[1]])))
    f <- function(row){
        if (all(is.na(row))) return(NA)
        if (all(!is.na(row)) && all(row == row[1])) return(row[1])
        row <- na.omit(row)
        if (!all(row == row[1])) return(".")
        return(NA)
    }
    result <- factor(apply(facMatrix, 1, f), levels = c(".", Levels))
    if (binary) result <- ifelse(result == ".", 0, 1)
    result
}
