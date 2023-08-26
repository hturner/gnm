# From issue #21

test_that("gnmFit handles linear fit provided optimal starting values", {
    Mean <- c(1237.1, 5605.55, 801.45, 2713.55, 570.4, 193.1, 97.2, 11.05,
              202.5, 7031.75, 2679.6, 1252.55, 735.6, 4088.05, 9818.7, 4486.45,
              3104.85, 1189.3, 217.6, 603.2, 28.45)
    VC <- c(33696.08, 296681.045, 24842.205, 31777.205, 1705.28, 950.48,
            2693.78, 244.205, 3026.42, 17578.125, 3.92, 8281.845, 1280.18,
            76479.605, 4665.78, 130101.005, 0.125, 9800, 18355.28, 1152,
            1618.805)
    dat <- data.frame(Mean=Mean,VC=VC)
    coeffs <- c(beta1 = 4792.94726157285, beta2 = 0.00366035757993686)
    
    # gnm with linear terms (calls glm.fit)
    fitglm <- gnm(VC ~ I(Mean^2) , family = Gamma(link = "identity"), 
                  data = data.frame(dat), start=coeffs, trace=TRUE)
    
    # gnm with equivalent "nonlin" term
    powfun <- function(x)
    {
        list(
            predictors=list(beta1 = 1, beta2 = 1),
            variables=list(substitute(x)),
            term=function(predictors,variables){
                paste( predictors[1],"+",predictors[2],"*",variables[1],"^2")
            }
        )
    }
    class(powfun) <- "nonlin"
    form <- VC ~ powfun(Mean)-1
    fitgnm <- gnm(formula = form, family = Gamma(link = "identity"), 
                  data = data.frame(dat), start=coeffs, trace=TRUE)

    expect_equal(coeffs, fitgnm$coefficients, ignore_attr = TRUE)
    # glm always does an extra iteration, even from previously converged fit
    expect_equal(coef(fitglm), coef(fitgnm), tolerance = 1e-6,
                 ignore_attr = TRUE)
})