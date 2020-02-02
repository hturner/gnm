context("implementation [Logistic]")

tol <- 1e-5

DNase1 <- subset(DNase, Run == 1)

## fit logistic model using nls
mod_nls <- nls(density ~ SSlogis( log(conc), Asym, xmid, scal ),
               data = DNase, subset = Run == 1)

## fit using basic nonlin terms
mod_basic <- gnm(density ~ -1 +
                     Mult(1, Inv(Const(1) + Exp(Mult(1 + offset(-log(conc)),
                                                     Inv(1))))),
                 start = c(NA, 0, 1), data = DNase1, verbose = FALSE)

## fit using Logistic()
Logistic <- function(x, inst = NULL){
    list(predictors = list(Asym = 1, xmid = 1, scal = 1),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             paste(predLabels[1], "/(1 + exp((", predLabels[2], "-",
                   varLabels[1], ")/", predLabels[3], "))")
         },
         start = function(theta){
             theta[3] <- 1
             theta
         }
    )
}
class(Logistic) <- "nonlin"
mod_logistic <- gnm(density ~ -1 + Logistic(log(conc)),
                    data = DNase1, verbose = FALSE)

test_that("logistic with gnm equivalent to nls", {
    expect_equivalent(unclass(coef(mod_basic)), coef(mod_nls), tol = tol)
    expect_equivalent(unclass(coef(mod_logistic)), coef(mod_nls), tol = tol)
    expect_equivalent(coef(mod_basic), coef(mod_logistic), tol = tol)
})
