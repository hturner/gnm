## as tolNorm2 method in rankMatrix but avoids validity checks
## much faster if need to do repeated rank calculations
quickRank <- function(X, tol = NULL) {
    sval <- svd(X, 0, 0)$d
    if (is.null(tol))
        sum(sval >= max(dim(X)) * .Machine$double.eps * sval[1])
    else
        sum(sval >= tol)
}
