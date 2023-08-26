tol <- 1e-4

test_that("Symm works with >2 factors", {
    # Table 8.5, Analysis of Ordinal Categorical Data, Agresti 2010 p242
    df8_5 <- data.frame(
        A = rep(1:3, each = 9),
        B = rep(rep(1:3, each = 3), times = 3),
        C = rep(1:3, times = 9),
        Freq = c(6, 4, 5, 3, 13, 10, 1, 8, 14,
                 2, 3, 2, 1, 3, 1, 2, 1, 2,
                 1, 0, 2, 0, 0, 0, 1, 1, 0)
    )
    ind <- expand.grid(1:3, 1:3, 1:3)
    ord_ind <- t(apply(ind, 1, sort))
    ref <- do.call("paste", c(as.data.frame(ord_ind), sep = ":"))
    expect_equal(with(df8_5, Symm(A, B, C)), as.factor(ref))
    
    # with named levels not in alphabetical order
    lev <- c("Red", "Blue", "Green")
    df8_5_nm <- data.frame(
        A = factor(rep(lev, each = 9), lev = lev),
        B = factor(rep(rep(lev, each = 3), times = 3), lev = lev),
        C = factor(rep(lev, times = 9), lev = lev),
        Freq = c(6, 4, 5, 3, 13, 10, 1, 8, 14,
                 2, 3, 2, 1, 3, 1, 2, 1, 2,
                 1, 0, 2, 0, 0, 0, 1, 1, 0)
    )
    res <- with(df8_5_nm, Symm(A, B, C))
    ord_lev <- apply(ord_ind, 2, function(x) lev[x])
    ref <- do.call("paste", c(as.data.frame(ord_lev), sep = ":"))
    ref <- factor(ref, unique(ref))
    expect_equal(with(df8_5_nm, Symm(A, B, C)), ref)
    
    # still works if factor elements are not in order
    scramble <- c(5, 12, 10, 23, 11, 15, 22, 6, 7, 21, 8, 19, 20,
                  25, 17, 16, 1, 18, 2, 14, 24, 26, 3, 13, 27, 4, 9)
    expect_equal(with(df8_5_nm[scramble,], Symm(A, B, C)), 
                 ref[scramble])
    
    # works for partial data - 3 x 3 x 2 array
    ref2 <- ref[df8_5_nm$C != "Green"]
    expect_equal(with(subset(df8_5_nm, C != "Green"), Symm(A, B, C)), 
                 droplevels(factor(ref2, levels(ref))))
    
    # works for partial data - 2 x 2 x 3 array
    ref3 <- ref[df8_5_nm$A != "Green" & df8_5_nm$B != "Green"]
    expect_equal(with(subset(df8_5_nm, A != "Green" & B != "Green"), 
                      Symm(A, B, C)), 
                 droplevels(factor(ref3, levels(ref))))
})
