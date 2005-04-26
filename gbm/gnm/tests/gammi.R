set.seed(1)
data(wheat)

bilinear3 <- gnm(yieldScaled ~ year + treatment +
                 Mult(year, treatment, multiplicity = 3),
                 family = gaussian, data = wheat)
bilinear3
