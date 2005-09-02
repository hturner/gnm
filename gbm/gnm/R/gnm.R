gnm <- function(formula, eliminate = NULL, constrain = NULL, family = gaussian,
                data = NULL, subset, weights, na.action,  method = "gnmFit",
                offset, start = NULL, control = gnmControl(...),
                verbose = TRUE, model = TRUE, x = FALSE, vcov = FALSE,
                termPredictors = FALSE, ...) {
    
    call <- match.call()
    
    modelTerms <- gnmTerms(formula, eliminate)
    falseFormula <- as.call(as.list(modelTerms))
    environment(falseFormula) <- environment(formula)
    
    modelData <- match.call(expand.dots = FALSE)
    argPos <- match(c("data", "subset", "weights", "na.action", "offset"),
                    names(modelData), 0)
    modelData <- as.call(c(as.name("model.frame"),
                           formula = falseFormula,
                           as.list(modelData)[argPos],
                           drop.unused.levels = TRUE))
    modelData <- eval(modelData, parent.frame())   
    attr(modelTerms, "variables") <- attr(attr(modelData, "terms"),
                                          "variables")

    if (!is.null(eliminate)) {
        if (!inherits(eliminate, "formula")) {
            stop("eliminate argument must be a formula")
        }
        elimTerms <- terms(eliminate)
        if (attr(elimTerms, "response") == 1) {
            stop("eliminate formula cannot have a response variable")
        }
        toElim <- attr(elimTerms, "factors")
        if (any(attr(attr(modelData, "terms"),
                     "dataClasses")[rownames(toElim)] != "factor"))
            stop("variables in 'eliminate' formula must be factors")
        elimCols <- model.matrix(update.formula(eliminate, ~ -1 + .),
                                                data = modelData)
        nElim <- ncol(elimCols)
        if (nrow(unique(elimCols)) > nElim)
            stop("'eliminate' formula is not equivalent to single factor")
    }
    else nElim <- 0
    
    if (method == "model.frame") {
        attr(modelData, "terms") <- attr(modelTerms, "terms")
        return(modelData)
    }
    else if (!method %in% c("gnmFit", "coefNames"))
        warning("method = ", method, " is not supported. Using \"gnmFit\".")
    
    y <- model.response(modelData, "numeric")
    nObs <- NROW(y)

    weights <- model.weights(modelData)
    if (!is.null(weights) & any(weights < 0))
        stop("negative weights are not allowed")
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

    if (is.empty.model(modelTerms)) {
        if (method == "coefNames") return(numeric(0))
        if (!family$valideta(offset))
            stop("invalid predictor values in empty model")
        mu <- family$linkinv(offset)
        if (!family$validmu(mu))
            stop("invalid fitted values in empty model")
        dmu <- family$mu.eta(offset)
        dev <- sum(family$dev.resids(y, mu, weights))
        modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs), mu,
                                                weights, dev))
        fit <- list(coefficients = numeric(0), eliminate = 0,
                    predictors = offset, fitted.values = mu, deviance = dev,
                    aic = modelAIC, iter = 0, conv = NULL,
                    weights = weights*dmu^2/family$variance(mu),
                    residuals = (y - mu)/dmu, df.residual = nObs, rank = 0)
        if (x) fit$x <- NULL
        if (vcov) fit$vcov <- NULL
        if (termPredictors) fit$termPredictors <- NULL
    }
    else {
        gnmEnvironment <- parent.frame()
        modelTools <- gnmTools(gnmEnvironment, modelTerms, modelData, x,
                               termPredictors)
        nParam <- length(modelTools$classID)

        if (method == "coefNames") return(names(modelTools$classID))

        if (is.null(constrain))
          constrain <- rep.int(FALSE, nParam)
        else {
            if (identical(constrain, "pick")) {
                if (!require(tcltk))
                    stop("constrain = \"pick\", and tcltk not installed")
                if (!require(relimp))
                    stop("the relimp package from CRAN needs to be installed")
                if (is.null(eliminate))
                    choice <- names(modelTools$classID)
                else
                    choice <- names(modelTools$classID)[-seq(nElim)]
                picked <- pickFrom(choice,
                              setlabels = "Coefficients to constrain",
                              title = "Constrain one or more gnm coefficients",
                              items.label = "Model coefficients:",
                              edit.setlabels = FALSE)
                constrain <- is.element(choice, picked)
                if (all(!constrain))
                    warning("no parameters were specified to constrain")
            }
            else if (is.numeric(constrain))
                constrain <- ifelse(is.element(seq(modelTools$classID) -
                                      nElim, constrain), TRUE, FALSE)
            else if (is.logical(constrain) & !is.null(eliminate))
                constrain <- c(rep(FALSE, nElim), constrain)
        }

        if (is.null(start))
            start <- rep.int(NA, nParam)
        else if (length(start) != nParam) {
            if (!is.null(eliminate) & length(start) == (nParam - nElim))
                start <- c(rep(NA, nElim), start)
            else
                stop("length(start) must either equal the no. of parameters\n",
                     "or the no. of non-eliminated parameters.")
        }
            
        
        fit <- gnmFit(modelTools, y, constrain, nElim, family, weights,
                       offset, nObs = nObs, start = start,
                       control = gnmControl(...), verbose, x, vcov,
                       termPredictors)
    }
    if (is.null(fit)) {
        warning("Algorithm failed - no model could be estimated", call. = FALSE)
        return()
    }
    fit <- c(list(call = call, formula = formula, constrain = constrain,
                  family = family, prior.weights = weights,
                  terms = attr(modelTerms, "terms"),
                  na.action = attr(modelData, "na.action"),
                  xlevels = .getXlevels(attr(modelData, "terms"), modelData),
                  y = y, offset = offset, control = control), fit)

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

