set.seed(1)
data(barley)

barleyRank2 <- gnm(y ~ -1 + Mult(site, variety, multiplicity = 2),
                    family = wedderburn, data = barley)
barleyRank2
