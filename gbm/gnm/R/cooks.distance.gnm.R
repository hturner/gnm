cooks.distance.gnm <- function(object, hat = hatvalues(object),
                                  dispersion = summary(object)$dispersion, ...){
    p <- object$rank
    res <- na.omit(residuals(object, type =
                             "pearson"))[object$prior.weights != 0]
    res <- naresid(object$na.action, res)
    (res/(1 - hat))^2 * hat/(dispersion * p)
}
