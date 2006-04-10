asGnm.lm <- function(object, ...) {
    lmExtra <- match(c("effects", "assign", "qr", "contrasts"), names(object))
    modelData <- model.frame(object)
    object[lmExtra] <- NULL
    object$call[[1]] <- as.name("gnm")
    constrain <- which(is.na(coef(object)))
    object <- c(list(formula = formula(object), eliminate = 0,
                     ofInterest = NULL, na.action = na.action(modelData),
                     constrain = constrain, constrainTo = numeric(constrain),
                     family = gaussian(), predictors = fitted.values(object),
                     deviance = deviance(object),
                     y = model.response(modelData)), object)
    object$weights <- object$prior.weights <- rep.int(1, length(object$y))
    object$aic <- 2 * object$rank +
        object$family$aic(object$y, object$weights, object$fitted.values,
                          object$weights, object$deviance)
    if (is.null(object$offset))
        object$offset <- rep.int(0, length(coef(object)))
    object$tolerance <- object$iterStart <- object$iterMax <- object$iter <-
        object$converged <- "Not available - model fitted by lm()"
    class(object) <- c("gnm", "glm", "lm")
    object
}
