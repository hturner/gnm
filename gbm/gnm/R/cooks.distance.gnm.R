cooks.distance.glm <- function(model, hat = hatvalues(model),
                                  dispersion = summary(model)$dispersion, ...){
    p <- model$rank
    res <- na.omit(residuals(model, type =
                             "pearson"))[model$prior.weights != 0]
    res <- naresid(model$na.action, res)
    (res/(1 - hat))^2 * hat/(dispersion * p)
}
