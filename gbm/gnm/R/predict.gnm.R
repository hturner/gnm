predict.gnm <- function (object, newdata = NULL,
                         type = c("link", "response", "terms"), se.fit = FALSE,
                         dispersion = NULL, terms = NULL, na.action = na.pass,
                         ...) {
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if(missing(newdata)) {
            pred <- switch(type, link = object$predictors,
                           response = object$fitted.values,
                           terms = object$termPredictors)
            if (!is.null(na.act))
                pred <- napredict(na.act, pred)
        }
        else {
            modelTerms <- delete.response(terms(object))
            modelData <- model.frame(modelTerms, newdata, na.action = na.action,
                                     xlev = object$xlevels)
            if (!is.null(offID <- attr(modelTerms, "offset")))
                offset <- eval(attr(modelTerms, "variables")[[offID + 1]],
                               newdata)
            else if (!is.null(object$offset))
                offset <- eval(object$call$offset, newdata)
            modelTools <- gnmTools(modelTerms, modelData, x = FALSE,
                                   termPredictors = type == "terms")
            varPredictors <-modelTools$varPredictors(coef(object))
            pred <- modelTools$predictor(varPredictors, term = type == "terms")
            if (!is.null(offset))  pred <- offset + pred
            switch(type, response = {pred <- family(object)$linkinv(pred)},
                   link = , terms = )
        }
    }
    else {
    }
    pred
}
