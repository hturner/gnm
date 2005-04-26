set.seed(1)
data(backPain)

oneDimensional <- update(noRelationship,
                         ~ . + Mult(pain - 1, x1 + x2 + x3 - 1),
                         startit = 3)
oneDimensional
