predict.gnm <- function (object, newdata = NULL,
                         type = c("link", "response", "terms"), se.fit = FALSE,
                         dispersion = NULL, terms = NULL, na.action = na.pass,
                         ...) {
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if(missing(newdata)) {
        pred <- switch(type, link = object$predictors,
                       response = object$fitted.values,
                       terms = {pred <- termPredictors(fit)
                                predc <- sweep(pred, 2, colMeans(pred))
                                const <- sum(pred[1,]) - sum(predc[1,])
                                structure(predc[, -1, drop = FALSE],
                                          constant = const)})
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
        modelTools <- gnmTools(modelTerms, modelData, x = FALSE)
        varPredictors <-modelTools$varPredictors(coef(object))
        pred <- modelTools$predictor(varPredictors, term = type == "terms")
        if (!is.null(offset))  pred <- offset + pred
        switch(type, response = {pred <- family(object)$linkinv(pred)},
               terms = {predc <- sweep(pred, 2, colMeans(termPredictors(fit)))
                        const <- sum(pred[1,]) - sum(predc[1,])
                        pred <- structure(predc[, -1], constant = const)},
               link = )
    }
    if (se.fit) {
        if (missing(newdata)) {
            V <- vcov(object, dispersion = dispersion)
            residual.scale <- as.vector(sqrt(attr(V, "dispersion")))
            X <- model.matrix(object)
            switch(type, link = {se.fit <- sqrt(diag(X %*% tcrossprod(V, X)))},
                   response = {
                       eta <- family(object)$linkfun(pred)
                       dX <- family(object)$mu.eta(eta) * X
                       se.fit <- sqrt(diag(dX %*% tcrossprod(V, dX)))},
                   terms = {
                       assign <- setdiff(unique(attr(X, "assign")), 0)
                       int <- attr(terms(object), "intercept")
                       X <- sweep(X, 2, colMeans(X))
                       se.fit <- matrix(, nr = nrow(pred), nc = ncol(pred))
                       for (i in assign + int)
                           se.fit[,i - int] <- sqrt(diag(X[, i] %*%
                                                   tcrossprod(V[i, i], X[, i])))
                   })
            if ( !is.null(na.act)) {
                se.fit <- napredict(na.act, se.fit)
            }
        }
        else {
        }
        dimnames(se.fit) <- dimnames(pred)
        pred <- list(fit = pred, se.fit = se.fit, residual.scale = residual.scale)
    }
    pred
}
