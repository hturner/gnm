library(gnm)
set.seed(1)
data(backPain)

library(nnet)
.incidence <- class.ind(backPain$pain)
.counts <- as.vector(t(.incidence))
.rowID <- factor(t(row(.incidence)))
backPain <- backPain[.rowID, ]
backPain$pain <- factor(rep(levels(backPain$pain), nrow(.incidence)),
                        levels = levels(backPain$pain), ordered = TRUE)

oneDimensional <- gnm(.counts ~ pain + Mult(pain - 1, x1 + x2 + x3 - 1),
                      eliminate = ~ .rowID, family = "poisson",
                      data = backPain, iterStart = 3)
oneDimensional
