context("datasets [voting]")

count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)

# standard Dref model
classMobility <- gnm(yvar ~ Dref(origin, destination),
                     constrain = "delta1",
                     family = binomial, data = voting, verbose = FALSE)

# separate weights for in and out of class 1
upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

socialMobility <- gnm(yvar ~ Dref(origin, destination,
                                  delta = ~ 1 + downward + upward),
                      constrain = "delta1",
                      family = binomial, data = voting, verbose = FALSE)

test_that("standard Dref model as expected for voting data", {
    expect_equal(round(deviance(classMobility), 2), 21.22)
    expect_equal(df.residual(classMobility), 19)
    p <- DrefWeights(classMobility)$origin["weight"]
    expect_equivalent(round(p, 2), 0.44)
})

test_that("modified Dref model as expected for voting data", {
    expect_equal(round(deviance(socialMobility), 2), 18.97)
    expect_equal(df.residual(socialMobility), 17)
    p <- DrefWeights(socialMobility)$origin[, "weight"]
    expect_equivalent(round(p, 2), c(0.40, 0.60, 0.39))
})
