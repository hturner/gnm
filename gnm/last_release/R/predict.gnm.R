#  Copyright (C) 2005, 2008, 2010, 2012, 2014 Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

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
        ## temporary fix; add new data to original data to get levels right
        modelData2 <- rbind(model.frame(object)[, colnames(modelData)],
                            modelData)
        if (length(offID <- attr(modelTerms, "offset")))
            offset <- eval(attr(modelTerms, "variables")[[offID + 1]],
                           newdata)
        else
            offset <- eval(object$call$offset, newdata)
        modelTools <- gnmTools(modelTerms, modelData2)
        varPredictors <- modelTools$varPredictors(parameters(object))
        pred <- modelTools$predictor(varPredictors, term = type == "terms")
        pred <- pred[-(1:length(object$y))]
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
        V <- vcov(object, dispersion = dispersion, with.eliminate = TRUE)
        residual.scale <- as.vector(sqrt(attr(V, "dispersion")))
        if (missing(newdata))
            X <- model.matrix(object)
        else
            X <- modelTools$localDesignFunction(parameters(object),
                                                varPredictors)
        covElim <- attr(V, "covElim")[object$eliminate, , drop = FALSE]
        switch(type,
               link = {
                   if (is.null(object$eliminate))
                       se.fit <- sqrt(diag(X %*% tcrossprod(V, X)))
                   else se.fit <-
                       sqrt(diag(X %*% tcrossprod(V, X)) +
                                2 * rowSums(X * covElim) +
                                    attr(V, "varElim")[object$eliminate])},
               response = {
                   eta <- na.omit(c(family(object)$linkfun(pred)))
                   d <- family(object)$mu.eta(eta)
                   if (is.null(object$eliminate))
                       se.fit <- sqrt(diag(X %*% tcrossprod(V, X)))
                   else se.fit <- sqrt(diag(X %*% tcrossprod(V, X)) +
                                       2*rowSums(X * covElim) +
                                       attr(V, "varElim")[object$eliminate])
                   se.fit <- d * se.fit},
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
