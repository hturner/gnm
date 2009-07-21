expandCategorical <- function(data, catvar, sep = ".", countvar = "count",
                              idvar = "id", sizevar = "size",
                              as.ordered = FALSE, group = TRUE) {

    cat <- interaction(data[catvar], sep = sep)
    ncat <- nlevels(cat)

    if (group == TRUE) {
        data <- data[, -match(catvar, names(data))]
        ord <- do.call("order", data)
        vars <- data[ord,]
        dupvars <- duplicated(vars)
        reps <- rle(dupvars)
        n <- sum(reps$values)
        ni <- reps$lengths[reps$values] + 1
        id <- factor(rep(seq(n), ni))

        counts <- as.data.frame(table(list(cat = cat[ord], id = id)))

        newData <- vars[which(!dupvars)[counts$id],]
        rownames(newData) <- NULL
        newData[c(catvar, idvar, countvar)] <- counts
        newData[[sizevar]] <- rep(ni, each = ncat)
    }
    else {
        n <- nrow(data)
        id <- gl(n, ncat)
        newData <- data[id, -match(catvar, names(data))]
        newData[[levels(interaction(catvar, sep = sep))]] <-
            gl(ncat, 1, n * ncat, labels = levels(cat), ordered = as.ordered)
        newData[[countvar]] <- as.vector(t(class.ind(cat)))
        newData[[idvar]] <- id
    }
    newData
}
