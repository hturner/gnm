predict.gnm <- function (object, newdata = NULL,
                         type = c("link", "response", "terms"), se.fit = FALSE,
                         dispersion = NULL, terms = NULL, na.action = na.pass,
                         ...) {
    type <- match.arg(type)
    na.act <- object$na.action
    if(missing(newdata)) {
        pred <- switch(type, link = object$predictors,
                       response = object$fitted.values,
                       terms = {pred <- termPredictors(object)
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
        else
            offset <- eval(object$call$offset, newdata)
        modelTools <- gnmTools(modelTerms, modelData)
        varPredictors <- modelTools$varPredictors(coef(object))
        pred <- modelTools$predictor(varPredictors, term = type == "terms")
        if (!is.null(offset))  pred <- offset + pred
        switch(type, response = {pred <- family(object)$linkinv(pred)},
               terms = {predc <- sweep(pred, 2, colMeans(termPredictors(object)))
                        const <- sum(pred[1,]) - sum(predc[1,])
                        pred <- structure(predc[, -1, drop = FALSE],
                                          constant = const)},
               link = )
    }
    if (se.fit) {
        V <- vcov(object, dispersion = dispersion)
        residual.scale <- as.vector(sqrt(attr(V, "dispersion")))
        if (missing(newdata))
            X <- model.matrix(object)
        else
            X <- modelTools$localDesignFunction(coef(object), varPredictors)
        switch(type, link = {se.fit <- sqrt(diag(X %*% tcrossprod(V, X)))},
               response = {
                   eta <- family(object)$linkfun(pred)
                   dX <- family(object)$mu.eta(eta) * X
                   se.fit <- sqrt(diag(dX %*% tcrossprod(V, dX)))},
               terms = {
                   assign <- setdiff(unique(modelTools$termAssign), 0)
                   int <- attr(terms(object), "intercept")
                   if (missing(newdata))
                       X <- sweep(X, 2, colMeans(X))
                   else
                       X <- sweep(X, 2, colMeans(model.matrix(object)))
                   se.fit <- matrix(, nr = nrow(pred), nc = ncol(pred))
                   for (i in assign + int)
                       se.fit[,i - int] <-
                           sqrt(diag(X[, i] %*% tcrossprod(V[i, i], X[, i])))
               })
        if ( !is.null(na.act)) {
            se.fit <- napredict(na.act, se.fit)
        }
        dimnames(se.fit) <- dimnames(pred)
        pred <- list(fit = pred, se.fit = se.fit,
                     residual.scale = residual.scale)
    }
    pred
}
