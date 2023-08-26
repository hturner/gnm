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
    expect_snapshot_value(doubleUnidiff, style = "serialize",
                          ignore_formula_env = TRUE)
    contr <- getContrasts(doubleUnidiff, 
                          rev(pickCoef(doubleUnidiff, ", class:vote")))
    expect_snapshot_value(contr, style = "serialize")
})
