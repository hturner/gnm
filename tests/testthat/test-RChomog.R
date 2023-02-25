context("implementation [multHomog]")

# Goodman, L. A. (1979) J. Am. Stat. Assoc., 74 (367), 537–552.

RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), family = poisson,
               data = occupationalStatus, verbose = FALSE)

test_that("RChomog model as expected for occupationalStatus data", {
    # Model (8) Table 7A
    pearson_chi_sq <- sum(na.omit(c(residuals(RChomog, type = "pearson")))^2)
    expect_equal(round(deviance(RChomog), 2), 32.56)
    expect_equal(round(pearson_chi_sq, 2), 31.21)
    expect_equal(df.residual(RChomog), 34)
})

# Chan, T.W. and Goldthorpe, J.H. (2004)  
# European Sociological Review, 20, 383–401.

###  Fit an association model with homogeneous row-column effects

### Set diagonal elements to NA (rather than fitting exactly)
dat <- as.data.frame(friend)
id <- with(dat, r == c)
dat[id,] <- NA
rc2 <- gnm(Freq ~ r + c + instances(MultHomog(r, c), 2),
           family = poisson, data = dat, iterStart = 0, verbose = FALSE)

test_that("RChomog2 model as expected for friend data", {
    # association models not reported in original paper
    pearson_chi_sq <- sum(na.omit(c(residuals(rc2, type = "pearson")))^2)
    expect_equal(round(deviance(rc2), 2), 1006.91)
    expect_equal(round(pearson_chi_sq, 2), 967.21)
    expect_equal(df.residual(rc2), 810)
})
