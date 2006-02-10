qrSolve2 <- function(A, b, rank = NULL, ...) {
    .C("dgelsy", m = as.integer(nrow(A)), n = as.integer(ncol(A)),
       nrhs = as.integer(ncol(b)), A, b,
       rcond = as.double(sqrt(.Machine$double.eps)))[[5]]
}
