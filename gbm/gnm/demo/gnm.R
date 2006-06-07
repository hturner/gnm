message("1. Load occupationalStatus data")
data(occupationalStatus)

message("2. Set seed as gnm returns random parameterization")
set.seed(1)

if (interactive()) {
    cat("\nType <Return> to start fitting models: ")
    readline()
}

message("3. Fit (linear) uniform association model, using Diag() to fit ",
        "diagonal effects")
Rscore <- scale(as.numeric(row(occupationalStatus)), scale = FALSE)
Cscore <- scale(as.numeric(col(occupationalStatus)), scale = FALSE)
Uniform <- gnm(Freq ~ origin + destination + Diag(origin, destination) + 
               Rscore:Cscore, family = poisson, data = occupationalStatus)
if (interactive()) {
    cat("\nType <Return> to print summary of Uniform: ")
    readline()
}
summary(Uniform)

message("4. Fit an association model using Mult() to fit separate row and ",
        "column effects")
RC <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
          Mult(origin, destination), family = poisson,
          data = occupationalStatus)
if (interactive()) {
    cat("\nType <Return> to print summary of RC: ")
    readline()
}
summary(RC)

message("5. Fit an association model using Nonlin() with MultHomog() plug-in\n",
        "to fit homogeneous row-column effects")
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               Nonlin(MultHomog(origin, destination)), family = poisson,
               data = occupationalStatus)
if (interactive()) {
    cat("\nType <Return> to print summary of RChomog: ")
    readline()
}
summary(RChomog)

message("6. Compare models using anova")
if (interactive()) {
    cat("\nType <Return> to print anova table : ")
    readline()
}
anova(Uniform, RChomog, RC)

message("7. Produce diagnostic plots for RChomog")
plot(RChomog)

message("8. Get simple constrasts of homogeneous row-column effects")
getContrasts(RChomog, grep("MultHomog", names(coef(RChomog))))

message("End of demo. \n",
 "See vignette("gnmOverview", package = \"gnm\") for full manual.")
