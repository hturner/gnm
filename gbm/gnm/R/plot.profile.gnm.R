plot.profile.gnm <- function (x, nseg, ...) 
{
    nulls <- sapply(x, is.null)
    if (all(nulls)) 
        return(NULL)
    x <- x[!nulls]
    pnames <- names(x)
    pnames <- pnames[!is.na(x[pnames])]
    nr <- ceiling(sqrt(length(pnames)))
    oldpar <- par(mfrow = c(nr, nr))
    on.exit(par(oldpar))
    for (nm in pnames) {
        z <- x[[nm]][[1]]
        parval <- x[[nm]][[2]][, nm]
        plot(parval, z, xlab = nm, ylab = "z", type = "n")
        if (sum(z == 0) == 1) 
            points(parval[z == 0], 0, pch = 3)
        splineVals <- spline(parval, z)
        lines(splineVals$x, splineVals$y)
    }
}
