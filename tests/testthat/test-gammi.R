tol <- 1e-4

# Vargas, M et al (2001). Interpreting treatment by environment interaction in 
# agronomy trials. Agronomy Journal 93, 949–960.

yield.scaled <- wheat$yield * sqrt(3/1000)
treatment <- interaction(wheat$tillage, wheat$summerCrop, wheat$manure,
                         wheat$N, sep = "")

mainEffects <- gnm(yield.scaled ~ year + treatment, data = wheat, 
                   verbose = FALSE)
svdStart <- residSVD(mainEffects, year, treatment, 3)
bilinear1 <- update(asGnm(mainEffects), . ~ . + 
                        Mult(year, treatment),
                    start = c(coef(mainEffects), svdStart[,1]))
bilinear2 <- update(bilinear1, . ~ . + 
                        Mult(year, treatment, inst = 2),
                    start = c(coef(bilinear1), svdStart[,2]))
bilinear3 <- update(bilinear2, . ~ . + 
                        Mult(year, treatment, inst = 3),
                    start = c(coef(bilinear2), svdStart[,3]))

test_that("bilinear model as expected for wheat data", {
    # check vs AMMI analysis of the T × E, end of Table 1
    txe <- anova(mainEffects, bilinear1, bilinear2, bilinear3)
    # year x treatment
    expect_equal(deviance(mainEffects), 279520, tolerance = tol)
    expect_equal(df.residual(mainEffects), 207)
    # diff for bilinear models
    expect_equal(txe$Deviance, c(NA, 151130, 39112, 36781), tolerance = tol)
    expect_equal(txe$Df, c(NA, 31, 29, 27))
    # "Deviations"
    expect_equal(deviance(bilinear3), 52497, tolerance = tol)
    expect_equal(df.residual(bilinear3), 120)
})
