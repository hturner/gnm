set.seed(1)
data(conformity)

A <- gnm(FCFF ~ Nonlin(Dref(FOPLF, MOPLM)) + FFCF + AGEM + FRMF +
         MRMM + MWORK - 1, family = gaussian, data = conformity)
A
F <- gnm(MCFM ~ Nonlin(Dref(MOPLM, FOPLF, formula = ~ 1 + MFCM)) + MFCM +
         AGEM + FRMF + MRMM + MWORK - 1, family = gaussian, data = conformity)
F
