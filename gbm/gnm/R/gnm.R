gnm <- function(formula, constrain = NULL, family = gaussian, data = NULL,
                subset, weights, na.action,  method = "gnm.fit", offset,
                start = NULL, control = gnm.control(...), model = TRUE,
                x = FALSE, vcov = FALSE, ...) {
    
    call <- match.call()
    
    modelTerms <- gnmTerms(formula)
    
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
    attr(modelTerms, "variables") <- attr(attr(modelData, "terms"),
                                          "variables")
    if (method == "model.frame") {
        attr(modelData, "terms") <- attr(modelTerms, "terms")
        return(modelData)
    }
    else if (!method %in% c("gnm.fit", "coef.names"))
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
                "multinomial response must be a factor")
        Y <- factor.incidence.matrix(y)
        resp.var.name <- names(modelData)[attr(attr(modelData, "terms"),
                                               "response")]
        if (is.null(colnames(Y))) colnames(Y) <- 1:ncol(Y)
        .rowID <- factor(t(row(Y)))
        assign(resp.var.name, factor(rep(colnames(Y), nrow(Y)),
                                     levels = levels(y),
                                     ordered = is.ordered(y)))
        .counts <- as.vector(t(Y))
        modelData <- modelData[.rowID, , drop = FALSE]
        modelData$.rowID <- .rowID
        modelData[[resp.var.name]] <- get(resp.var.name)
        newCall <- call
        newCall$formula <- update.formula(formula,
                                          .counts ~ -1 + .rowID + .)
        newCall$family <- as.name("poisson")
        newCall$data <- as.name("modelData")
        if (!is.null(call$weights))
          newCall$weights <- rep(weights * rowSums(Y), rep(ncol(Y), nrow(Y)))
        if (!is.null(constrain)) {
            if (constrain == "pick") {
                if (!require(tcltk))
                    stop("constrain = \"pick\", and tcltk not installed")
                if (!require(relimp))
                    stop("the relimp package from CRAN needs to be installed")
                constrain <- is.element(names(modelTools$classIndex),
                                        pickFrom(names(modelTools$classIndex)
                                                 [-seq(nrow(Y))]))
                if (all(!constrain))
                    warning("no parameters were specified to constrain")
            }
            if (is.numeric(constrain))
                newCall$constrain <- constrain + nrow(Y)
            else
                newCall$constrain <- c(rep(FALSE, nrow(Y)), constrain)
        }
        newCall$subset <- NULL
        result <- eval(newCall)
        result$original.call <- call
        result$auxiliary[seq(nrow(Y))] <- TRUE
        return(result)
}

    if (is.empty.model(modelTerms)) {
        if (method == "coef.names") return(numeric(0))
        if (!family$valideta(offset))
            stop("Invalid predictor values in empty model")
        mu <- family$linkinv(offset)
        if (!family$validmu(mu))
            stop("Invalid fitted values in empty model")
        dmu <- family$mu.eta(offset)
        dev <- sum(family$dev.resids(y, mu, weights))
        modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs), mu,
                                                weights, dev))
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
        modelTools <- gnmTools(modelTerms, modelData, x, family, weights,
                               offset)

        if (method == "coef.names") return(names(modelTools$classIndex))

        if (is.null(constrain))
          constrain <- rep.int(FALSE, length(modelTools$classIndex))
        else {
            if (constrain == "pick") {
                if (!require(tcltk))
                    stop("constrain = \"pick\", and tcltk not installed")
                if (!require(relimp))
                    stop("the relimp package from CRAN needs to be installed")
                constrain <- is.element(names(modelTools$classIndex),
                                        pickFrom(names(modelTools$classIndex)))
                if (all(!constrain))
                    warning("no parameters were specified to constrain")
            }
            if (is.numeric(constrain))
                constrain <- ifelse(is.element(seq(modelTools$classIndex),
                                               constrain), TRUE, FALSE)
        }

        if (is.null(start)) start <- rep.int(NA, length(modelTools$classIndex))
        
        fit <- gnm.fit(modelTools, y, constrain, family, weights,
                       offset, nObs = nObs, start = start,
                       control = gnm.control(...), x, vcov)
    }
    fit <- c(list(call = call, formula = formula, constrain = constrain,
                  family = family, prior.weights = weights,
                  terms = attr(modelTerms, "terms"),
                  na.action = attr(modelData, "na.action"),
                  xlevels = .getXlevels(attr(modelData, "terms"), modelData),
                  y = y, offset = offset, control = control), fit,
             list(auxiliary = rep(FALSE, length(coef(fit)))))

    asY <- c("predictors", "fitted.values", "residuals", "prior.weights",
             "weights", "y", "offset")
    if (inherits(data, "table") & !is.empty.model(modelTerms)) {
        fit[asY] <- lapply(fit[asY], function(x) {
            as.table(tapply(x,
                            as.data.frame(data)[names(dimnames(data))], sum))})
    }
    else
        fit[asY] <- lapply(fit[asY],
                           function(x, y) structure(x, names = names(y)), y)
    
    if (model) {
        attr(modelData, "terms") <- attr(modelTerms, "terms")
        fit$model <- modelData
    }
    class(fit) <- c("gnm", "glm", "lm")
    fit
}

