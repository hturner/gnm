context("model spec [Symm]")

tol <- 1e-4

test_that("confint works within function call", {
    # https://github.com/hturner/gnm/issues/11
    f <- function(d) {
        fit <- gnm(Freq ~ vote, family = poisson, data = d)
        confint(fit)
    }
    expect_known_value(f(as.data.frame(cautres)),
                       file = test_path("outputs/confint.rds"))
})
