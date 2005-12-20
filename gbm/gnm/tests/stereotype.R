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

noRelationship <- gnm(.counts ~ pain, eliminate = ~ .rowID,
                      family = "poisson", data = backPain)

oneDimensional <- update(noRelationship,
                         ~ . + Mult(pain - 1, x1 + x2 + x3 - 1))
oneDimensional
anova(oneDimensional)
anova(noRelationship, oneDimensional)
coef(oneDimensional)
cooks.distance(oneDimensional)
labels(oneDimensional)
getContrasts(oneDimensional, 11:6)
model.frame(oneDimensional)
model.matrix(oneDimensional)
plot(oneDimensional, ask = FALSE)
residuals(oneDimensional)
rstandard(oneDimensional)
summary(oneDimensional)
termPredictors(oneDimensional)
variable.names(oneDimensional)
vcov(oneDimensional)

