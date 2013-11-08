#  Copyright (C) 2006, 2009, 2013 Heather Turner
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

expandCategorical <- function(data, catvar, sep = ".", countvar = "count",
                              idvar = "id", as.ordered = FALSE, group = TRUE) {

    cat <- interaction(data[catvar], sep = sep)
    ncat <- nlevels(cat)

    if (group == TRUE) {
        data <- data[, -match(catvar, names(data))]
        ord <- do.call("order", data)
        vars <- data[ord,]
        dupvars <- duplicated(vars)
        reps <- rle(dupvars)
        n <- sum(!dupvars)
        ni <- rep(1, n)
        grp.ni <- reps$lengths[reps$values] + 1
        grps <- cumsum(reps$lengths[!reps$values])[seq(grp.ni)]
        ni[grps] <- grp.ni
        id <- factor(rep(seq(n), ni))

        counts <- as.data.frame(table(list(cat = cat[ord], id = id)))

        newData <- vars[which(!dupvars)[counts$id],]
        rownames(newData) <- NULL
        newData[c(catvar, idvar, countvar)] <- counts
    }
    else {
        n <- nrow(data)
        id <- gl(n, ncat)
        newData <- data[id, -match(catvar, names(data)), drop = FALSE]
        newData[[levels(interaction(catvar, sep = sep))]] <-
            gl(ncat, 1, n * ncat, labels = levels(cat), ordered = as.ordered)
        newData[[countvar]] <- as.vector(t(class.ind(cat)))
        newData[[idvar]] <- id
    }
    newData
}
