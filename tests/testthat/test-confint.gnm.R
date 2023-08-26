tol <- 1e-4

test_that("confint works within function call", {
    # https://github.com/hturner/gnm/issues/11
    f <- function(d) {
        fit <- gnm(Freq ~ vote, family = poisson, data = d)
        confint(fit)
    }
    expect_snapshot_value(f(as.data.frame(cautres)), style = "json2")
})
