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
        newData <- data[id, -match(catvar, names(data))]
        newData[[levels(interaction(catvar, sep = sep))]] <-
            gl(ncat, 1, n * ncat, labels = levels(cat), ordered = as.ordered)
        newData[[countvar]] <- as.vector(t(class.ind(cat)))
        newData[[idvar]] <- id
    }
    newData
}
