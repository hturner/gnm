gnm <- function(formula, constrain = NULL, family = gaussian, data = NULL,
                subset, weights, na.action,  method = "gnm.fit", offset,
                start = NULL, control = gnm.control(...), model = TRUE,
                x = FALSE, vcov = FALSE, ...) {
    
    call <- match.call()
    
    modelTerms <- gnmTerms(formula)
    
    if (inherits(data, "table")) {
        tableFormula <- paste(names(dimnames(data)), collapse = "+")
        data <- as.data.frame(data)
    }
    else tableFormula <- NULL
    
    modelData <- match.call(expand.dots = FALSE)
    argPos <- match(c("data", "subset", "weights", "na.action", "offset"),
                    names(modelData), 0)
    modelData <- as.call(c(as.name("model.frame"),
                           formula = as.call(as.list(modelTerms)),
                           as.list(modelData)[argPos],
                           drop.unused.levels = TRUE))
    if (!is.null(attr(modelTerms, "extraData")))
        modelData <- eval(modelData, attr(modelTerms, "extraData"),
                          parent.frame())
    else
        modelData <- eval(modelData, parent.frame())   
    attr(modelTerms, "variables") <- attr(attr(modelData, "terms"), "variables")
    if (method == "model.frame") {
        attr(modelData, "terms") <- attr(modelTerms, "terms")
        return(modelData)
    }
    else if (method != "gnm.fit")
        warning("method = ", method, " is not supported. Using \"gnm.fit\".")
    
    y <- model.response(modelData, "numeric")
    nObs <- NROW(y)

    weights <- model.weights(modelData)
    if (!is.null(weights) & any(weights < 0))
        stop("Negative weights are not allowed")
    if (is.null(weights))
        weights <- rep.int(1, nObs)
    offset <- model.offset(modelData)
    if (is.null(offset))
        offset <- rep.int(0, nObs)

    if (!missing(family))
        if (!exists(as.character(call$family)))
            stop("family object \"", as.character(call$family), "\" not found")
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
    
    if (family$family == "binomial") {
        if (NCOL(y) == 1 & is.factor(y))
                y <- y != levels(y)[1]
        else if (NCOL(y) == 2) {
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            weights <- weights * n
        }
    }

    if (family$family == "multinomial") {
        if (!is.null(call$offset)) stop(
                "offset cannot be used with multinomial logit models") 
        if (!is.factor(y)) stop(
                "multinomial response data must be factor or matrix")
        y <- factor.incidence.matrix(y)
        resp.var.name <- names(modelData)[attr(attr(modelData, "terms"),
                                               "response")]
        if (is.null(colnames(y))) colnames(y) <- 1:ncol(y)
        .rowID <- factor(t(row(y)))
        assign(resp.var.name, factor(rep(colnames(y), nrow(y))))
        .counts <- as.vector(t(y))
        modelData <- modelData[.rowID, , drop = FALSE]
        modelData$.rowID <- .rowID
        modelData[[resp.var.name]] <- get(resp.var.name)
        newCall <- call
        newCall$formula <- update.formula(formula,
                                          .counts ~ -1 + .rowID + .)
        newCall$family <- as.name("poisson")
        newCall$data <- as.name("modelData")
        if (!is.null(call$weights))
          newCall$weights <- rep(weights * rowSums(y), rep(ncol(y), nrow(y)))
        if (!is.null(constrain)) {
          if (is.list(constrain))
            newCall$constrain <- lapply(constrain, "+", nrow(y))
          else
            newCall$constrain <- c(rep(FALSE, nrow(y)), constrain)
        }
        newCall$subset <- NULL
        result <- eval(newCall)
        result$original.call <- call
        result$auxiliary[seq(nrow(y))] <- TRUE
        return(result)
}

    if (is.empty.model(modelTerms)) {
        if (!valideta(offset))
            stop("Invalid predictor values in empty model")
        mu <- family$linkinv(offset)
        if (!validmu(mu))
            stop("Invalid fitted values in empty model")
        dmu <- family$mu.eta(eta)
        dev <- sum(family$dev.resids(y, mu, weights))
        modelAIC <- family$aic(y, rep.int(1, nObs), mu, weights, dev) +
        2 * attr(VCOV, "rank")
        fit <- list(coefficients = numeric(0), predictors = offset,
                    fitted.values = mu, deviance = dev, aic = modelAIC,
                    iter = 0, conv = NULL,
                    weights = weights*dmu^2/family$variance(mu),
                    residuals = (y - mu)/dmu,
                    df.residual = nObs, rank = 0)
        if(x) fit$x <- NULL
        if (vcov) fit$vcov <- NULL
    }
    else {
        modelTools <- gnmTools(modelTerms, modelData, x)

        if (is.null(constrain))
          constrain <- rep.int(FALSE, length(modelTools$factorAssign))
        else if (is.list(constrain))
          constrain <- ifelse(is.element(seq(modelTools$factorAssign),
                                         constrain), TRUE, FALSE)

        fit <- gnm.fit(modelTools, y, constrain, family, weights,
                        offset, nObs = nObs, start = start,
                        control = gnm.control(...), x, vcov)
    }
    fit <- c(list(call = call, formula = formula, constrain = constrain,
                  family = family, prior.weights = weights,
                  terms = attr(modelTerms, "terms"),
                  na.action = attr(modelData, "na.action"),
                  xlevels = .getXlevels(modelTerms, modelData), y = y,
                  offset = offset, control = control), fit,
             list(auxiliary = rep(FALSE, length(coef(fit)))))
    
    if (!is.null(tableFormula)) {
        toTable <- c("predictors", "fitted.values", "residuals",
                     "prior.weights", "weights", "y", "offset")
        fit[toTable] <- lapply(toTable, function(x)
                               xtabs(as.formula(paste(x, "~", tableFormula)),
                                     data = cbind(modelData, fit[toTable])))
    }
    if (model) {
        attr(modelData, "terms") <- attr(modelTerms, "terms")
        fit$model <- modelData
    }
    class(fit) <- c("gnm", "glm", "lm")
    fit
}

