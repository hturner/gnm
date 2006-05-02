library(gnm)
set.seed(1)
data(barley)

biplotModel <- gnm(y ~ -1 + Mult(site, variety, multiplicity = 2),
                   family = wedderburn, data = barley)
biplotModel
