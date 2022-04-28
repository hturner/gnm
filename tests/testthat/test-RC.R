context("datasets [mentalHealth]")

# set seed to fix sign of coef
suppressWarnings(RNGversion("3.0.0")) 
set.seed(1)

# Agresti A (2002).Categorical Data Analysis. 2nd edition

mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)

RC1model <- gnm(count ~ SES + MHS +
                    Mult(-1 + SES, -1 + MHS),
                family = poisson, data = mentalHealth, verbose = FALSE)

test_that("RC model as expected for mentalHealth data", {
    # compare vs results in sec 9.6.2
    pearson_chi_sq <- sum(na.omit(c(residuals(RC1model, type = "pearson")))^2)
    expect_equal(round(pearson_chi_sq, 1), 3.6)
    expect_equal(df.residual(RC1model), 8)
    # normalize as in Agresti's eqn 9.15
    rowProbs <- with(mentalHealth, tapply(count, SES, sum) / sum(count))
    colProbs <- with(mentalHealth, tapply(count, MHS, sum) / sum(count))
    mu <- getContrasts(RC1model, pickCoef(RC1model, "[.]SES"),
                       ref = rowProbs, scaleRef = rowProbs,
                       scaleWeights = rowProbs)
    nu <- getContrasts(RC1model, pickCoef(RC1model, "[.]MHS"),
                       ref = colProbs, scaleRef = colProbs,
                       scaleWeights = colProbs)
    # change of scale
    expect_equal(round(-mu$qvframe$Estimate, 2),  
                 c(-1.11, -1.12, -0.37, 0.03, 1.01, 1.82))
    expect_equal(round(-nu$qvframe$Estimate, 2),  
                 c(-1.68, -0.14, 0.14, 1.41))
    # association parameter
    rowScores <- coef(RC1model)[10:15]
    colScores <- coef(RC1model)[16:19]
    rowScores <- rowScores - sum(rowScores * rowProbs)
    colScores <- colScores - sum(colScores * colProbs)
    beta1 <- sqrt(sum(rowScores^2 * rowProbs))
    beta2 <- sqrt(sum(colScores^2 * colProbs))
    expect_equal(round(beta1 * beta2, 2), 0.17)
})
