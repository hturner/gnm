# lwork = max(min(m, n) + 3 * n + 1, 2 * min(m, n) + nrhs)
qrSolve2 <- function(A, b, rank = NULL, ...) {
    m <- as.integer(NROW(A))
    n <- as.integer(NCOL(A))
    bMissing <- FALSE
    if (missing(b)){
        if (m != n) stop("only square matrices can be inverted")
        b <- diag(m)
        bMissing <- TRUE
    }
    nrhs <- as.integer(NCOL(b))
    ans <- matrix(nrow = n, ncol = nrhs)
    if (m < 1)
        stop("A is empty")
    if (NROW(b) != m)
        stop("NROW(b) != NROW(A)")
    if (m < n)
        b <- rbind(as.matrix(b), matrix(nrow = n - m, ncol = nrhs))
    storage.mode(A) <- storage.mode(b) <- storage.mode(ans) <- "double"
    work <- double(1)
    lwork <- as.integer(-1)
    result <- .C("dgelsy", m = m, n = n, nrhs = nrhs, A = A, b = b,
                 rcond = as.double(3e-17), rank = integer(1), work = work,
                 lwork = lwork, ans = ans, NAOK = TRUE)
    lwork <- as.integer(result$work[1])
    work <- double(lwork)
    result <- .C("dgelsy", m = m, n = n, nrhs = nrhs, A = A, b = b,
                 rcond = as.double(100 * .Machine$double.eps),
                 rank = integer(1), work = work, lwork = lwork, ans = ans,
                 NAOK = TRUE)
    result <- structure(drop(result$ans), rank = result$rank)
    if (is.matrix(result)) {
        rownames(result) <- colnames(A)
        colnames(result) <- if (bMissing) rownames(A)
        else if (is.matrix(b)) colnames(b)
        else names(b)
    } else names(result) <- colnames(A)
    result
}
