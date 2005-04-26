set.seed(1)
data(mentalHealth)

mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)
RC1model <- gnm(count ~ SES + MHS +
                Mult(-1 + SES, -1 + MHS),
                family = poisson, data = mentalHealth)
RC1model
