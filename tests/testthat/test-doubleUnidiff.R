context("datasets [cautres]")

# set seed to compare to saved values (not all identifiable)
suppressWarnings(RNGversion("3.0.0")) 
set.seed(1)

test_that("double unidiff model as expected for cautres data", {
    doubleUnidiff <- gnm(Freq ~ election*vote + election*class*religion +
                             Mult(Exp(election), religion:vote) +
                             Mult(Exp(election), class:vote),
                         family = poisson, data = cautres, verbose = FALSE)
    doubleUnidiff$family$dispersion <- 1
    expect_equal(round(deviance(doubleUnidiff), 2), 133.04)
    expect_known_value(doubleUnidiff,
                       file = test_path("outputs/doubleUnidiff.rds"))
    contr <- getContrasts(doubleUnidiff, 
                          rev(pickCoef(doubleUnidiff, ", class:vote")))
    expect_known_value(contr,
                       file = test_path("outputs/doubleUnidiff-contrasts.rds"))
})
