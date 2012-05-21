predict.gnm <- function (object, newdata = NULL,
                         type = c("link", "response", "terms"), se.fit = FALSE,
                         dispersion = NULL, terms = NULL,
                         na.action = na.exclude, ...) {
    type <- match.arg(type)
    if (type == "terms") {
        hasintercept <- attr(object$terms, "intercept") > 0L
        if (is.null(terms))
            terms <- c("(eliminate)"[!is.null(object$eliminate)],
                       attr(object$terms, "term.labels"))
    }
    if(missing(newdata)) {
        pred <- switch(type, link = object$predictors,
                       response = object$fitted.values,
                       terms = {pred <- termPredictors(object)
                                ## see 6.3.6 white book & predict.lm
                                if (hasintercept) {
                                    predc <- sweep(pred, 2, colMeans(pred))
                                    const <- sum(pred[1,]) - sum(predc[1,])
                                    structure(predc[, terms, drop = FALSE],
                                              constant = const)
                                } else structure(pred[, terms, drop = FALSE],
                                                 constant = 0)})
        if (!is.null(na.act <- object$na.action)){
            pred <- napredict(na.act, pred)
        }
        if (!inherits(pred, "matrix") && !is.null(object$table.attr))
            attributes(pred) <- object$table.attr
    }
    else {
        modelTerms <- delete.response(terms(object))
        modelData <- model.frame(modelTerms, newdata, na.action = na.action,
                                 xlev = object$xlevels)
        if (length(offID <- attr(modelTerms, "offset")))
            offset <- eval(attr(modelTerms, "variables")[[offID + 1]],
                           newdata)
        else
            offset <- eval(object$call$offset, newdata)
        modelTools <- gnmTools(modelTerms, modelData)
        varPredictors <- modelTools$varPredictors(parameters(object))
        pred <- modelTools$predictor(varPredictors, term = type == "terms")
        names(pred) <- rownames(modelData)
        if (!is.null(offset))  pred <- offset + pred
        switch(type, response = {pred <- family(object)$linkinv(pred)},
               terms = {predc <- sweep(pred, 2, colMeans(termPredictors(object)))
                        const <- sum(pred[1,]) - sum(predc[1,])
                        pred <- structure(predc[, terms, drop = FALSE],
                                              constant = const)},
               link = )
        if (!is.null(na.act <- attr(modelData, "na.action")))
            pred <- napredict(na.act, pred)
    }
    if (se.fit) {
        V <- vcov(object, dispersion = dispersion)
        residual.scale <- as.vector(sqrt(attr(V, "dispersion")))
        if (missing(newdata))
            X <- model.matrix(object)
        else
            X <- modelTools$localDesignFunction(parameters(object),
                                                varPredictors)
        e <- object$eliminate
        switch(type,
               link = {
                   if (is.null(e))
                       se.fit <- sqrt(diag(X %*% tcrossprod(V, X)))
                   else se.fit <- sqrt(diag(X %*% tcrossprod(V, X)) +
                                       2 * rowSums(X * attr(V, "covElim")[e,]) +
                                       attr(var, "varElim")[e])},
               response = {
                   eta <- na.omit(c(family(object)$linkfun(pred)))
                   d <- family(object)$mu.eta(eta)
                   dX <- d * X
                   if (is.null(e))
                       se.fit <- sqrt(diag(dX %*% tcrossprod(V, dX)))
                   else se.fit <- sqrt(diag(dX %*% tcrossprod(V, dX)) +
                                       2*rowSums(dX * attr(V, "covElim")[e,]) +
                                       d * attr(var, "varElim")[e])},
               terms = {
                   if (missing(newdata)) {
                       assign <- split(seq(ncol(X)), attr(X, "assign"))
                       X <- sweep(X, 2, colMeans(X))
                   }
                   else {
                       M <- model.matrix(object)
                       assign <- split(seq(ncol(X)), attr(M, "assign"))
                       X <- sweep(X, 2, colMeans(M))
                   }
                   se.fit <- matrix(, nrow = nrow(X), ncol = length(terms))
                   s <- 0
                   for (i in match(terms, colnames(pred))) {
                       s <- s + 1
                       t <- assign[[i]]
                       se.fit[, s] <-
                           sqrt(diag(X[, t] %*% tcrossprod(V[t, t], X[, t])))
                   }
               })
        if (!is.null(na.act)) {
            se.fit <- napredict(na.act, se.fit)
        }
        if (inherits(pred, "table"))
            attributes(se.fit) <- object$table.attr
        else
            attributes(se.fit) <- attributes(pred)
        pred <- list(fit = pred, se.fit = se.fit,
                     residual.scale = residual.scale)
    }
    pred
}
