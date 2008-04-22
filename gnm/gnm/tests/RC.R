library(gnm)
set.seed(1)
data(mentalHealth)

mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)
RC1model <- gnm(count ~ SES + MHS +
                Mult(-1 + SES, -1 + MHS),
                family = poisson, data = mentalHealth)

print(RC1model$deviance, digits = 10)
print(RC1model$df)

RC1modelOffset <- gnm(count ~ SES + MHS + Const(3) + 
                      Mult(-1 + SES, -1 + MHS),
                      family = poisson, data = mentalHealth)
coef(RC1model)[1] - coef(RC1modelOffset)[1]
