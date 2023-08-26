tol <- 1e-4

library(MASS)
example(birthwt, echo = FALSE)
library(nnet)
bwt.mu <- multinom(low ~ ., data = bwt, trace = FALSE)

## Equivalent using gnm - include unestimable main effects in model so
## that interactions with low0 automatically set to zero, else could use
## 'constrain' argument.
bwtLong <- expandCategorical(bwt, "low", group = FALSE)
bwt.po <- gnm(count ~  low*(. - id), eliminate = id, family = "poisson", 
              data = bwtLong, verbose = FALSE)

## Equivalent using glm
bwt.po2 <- glm(formula = count ~ -1 + id + low * (. -id), family = "poisson", 
               data = bwtLong)

test_that("gnm agrees with multinom", {
    cf0 <- coef(bwt.mu)
    cf1 <- na.omit(coef(bwt.po))
    expect_equal(cf0, cf1, tolerance = tol, ignore_attr = TRUE)
    expect_equal(deviance(bwt.mu), deviance(bwt.po), tolerance = tol,
                 ignore_attr = TRUE)
})

test_that("gnm agrees with glm", {
    cf1 <- coef(bwt.po)
    all_coef1 <- c(attr(cf1, "eliminated"), cf1)
    expect_equal(all_coef1, coef(bwt.po2), tolerance = tol, ignore_attr = TRUE)
})


