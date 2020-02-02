context("datasets [barley]")

# set seed to fix sign
suppressWarnings(RNGversion("3.0.0")) 
set.seed(1)

# Gabriel, K R (1998). Generalised bilinear regression. Biometrika 85, 689–700.

test_that("biplot model as expected for barley data", {
    biplotModel <- gnm(y ~ -1 + instances(Mult(site, variety), 2),
                       family = wedderburn, data = barley, verbose = FALSE)
    expect_known_value(biplotModel,
                       file = test_path("outputs/biplotModel.rds"))
    # rotate and scale fitted predictors
    barleyMatrix <- xtabs(biplotModel$predictors ~ site + variety,
                          data = barley)
    barleySVD <- svd(barleyMatrix)
    A <- sweep(barleySVD$u, 2, sqrt(barleySVD$d), "*")[, 1:2]
    B <- sweep(barleySVD$v, 2, sqrt(barleySVD$d), "*")[, 1:2]
    rownames(A) <- levels(barley$site)
    rownames(B) <- levels(barley$variety)
    colnames(A) <- colnames(B) <- paste("Component", 1:2)
    # compare vs matrices in Gabriel (1998): allow for sign change
    # 3rd element in fit is 1.425 vs 1.42 in paper
    expect_equivalent(round(A, 2), 
                      matrix(c(4.19, 2.76, 1.43, 1.85, 1.27, 
                               1.16, 1.02, 0.65, -0.15, 
                               -0.39, -0.34, -0.05, 0.33, 0.16, 
                               0.4, 0.73, 1.46, 2.13), nrow = 9))
    expect_equivalent(round(B, 2), 
                      matrix(c(-2.07, -3.06, -2.96, -1.81, -1.56,
                               -1.89, -1.18, -0.85, -0.97, -0.60,
                               -0.97, -0.51, -0.33, -0.50, -0.08, 
                               1.08, 0.41, 1.15, 1.27, 1.40), nrow = 10))
    # chi-square statistic approx equal to that reported
    expect_equal(round(sum(residuals(biplotModel, type = "pearson")^2)), 54)
    expect_equal(df.residual(biplotModel), 56)
})

