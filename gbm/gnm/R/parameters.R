parameters <- function(object) {
    coefs <- coef(object)
    coefs[object$constrain[,1]] <- object$constrain[,2]
    coefs
}
