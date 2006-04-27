library(gnm)
set.seed(1)
data(barley)

biplotModel <- gnm(y ~ -1 + instances(Mult(site, variety), 2),
                   family = wedderburn, data = barley)
biplotModel
