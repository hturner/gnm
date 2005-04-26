set.seed(1)
data(occupational.status)

RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               Nonlin(MultHomog(origin, destination)), family = poisson,
               data = occupational.status)
RChomog
