asGnm.glm <- function(object, ...) {
    glmExtra <- match(c("effects", "R", "qr", "null.deviance", "df.null",
                        "boundary", "control", "contrasts"), names(object))
    modelData <- model.frame(object)
    object[glmExtra] <- NULL
    object$call[[1]] <- as.name("gnm")
    constrain <- which(is.na(coef(object)))
    object$terms <- gnmTerms(object$formula, data = modelData)
    object <- c(list(eliminate = 0, ofInterest = NULL,
                     na.action = na.action(modelData), constrain = constrain,
                     constrainTo = numeric(length(constrain))),
                object)
    names(object)[match("linear.predictors", names(object))] <- "predictors"
    if (is.null(object$offset))
        object$offset <- rep.int(0, length(coef(object)))
    object$tolerance <- object$iterStart <- object$iterMax <-
        "Not available - model fitted by glm()"
    class(object) <- c("gnm", "glm", "lm")
    object
}
