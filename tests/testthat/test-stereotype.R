context("datasets [backPain]")

tol <- 1e-4

backPainLong <- expandCategorical(backPain, "pain")

## stereotype model
stereotype <- gnm(count ~ pain + Mult(pain, x1 + x2 + x3), 
                  eliminate = id, family = "poisson",
                  data = backPainLong, verbose = FALSE)

test_that("sterotype model as expected for backPain data", {
    # Obtain number of parameters and log-likelihoods for equivalent
    # "Six groups: one-dimensional" multinomial model presented in Table 5
    # maximised log-likelihood
    size <- tapply(backPainLong$count, backPainLong$id, sum)[backPainLong$id]
    expect_equal(round(sum(stereotype$y * log(stereotype$fitted/size)), 2),
                 -151.55)
    # number of parameters
    expect_equivalent(stereotype$rank - nlevels(stereotype$eliminate), 12)
})
