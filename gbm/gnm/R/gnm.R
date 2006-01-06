gnm <- function(formula, eliminate = NULL, constrain = NULL, family = gaussian,
                data = NULL, subset, weights, na.action,  method = "gnmFit",
                offset, start = NULL, control = gnmControl(...),
                verbose = TRUE, model = TRUE, x = FALSE, vcov = FALSE,
                termPredictors = FALSE, ...) {
    
    call <- match.call()
    
    modelTerms <- gnmTerms(formula, substitute(eliminate))
    modelData <- match.call(expand.dots = FALSE)
    argPos <- match(c("data", "subset", "weights", "na.action", "offset"),
                    names(modelData), 0)
    modelData <- as.call(c(as.name("model.frame"),
                           formula = modelTerms,
                           as.list(modelData)[argPos],
                           drop.unused.levels = TRUE))
    modelData <- eval(modelData, parent.frame())
    
    if (method == "model.frame") {
        attr(modelData, "terms") <- attr(modelTerms, "terms")
        return(modelData)
    }
    else if (!method %in% c("gnmFit", "coefNames", "model.matrix")) {
        warning("method = ", method, " is not supported. Using \"gnmFit\".",
                call. = FALSE)
        method <- "gnmFit"
    }

    if (!missing(eliminate)) {
        Elim <- modelData[[attr(attr(modelTerms, "terms"), "term.labels")[1]]]
        if (!is.factor(Elim))
            stop("variables in 'eliminate' formula must be factors")
        nElim <- nlevels(Elim)
    }
    else nElim <- 0
    
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
        else if (method == "model.matrix")
            return(model.matrix(formula, data = modelData))
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
        if (x) fit <- c(fit, x = NULL)
        if (vcov) fit <- c(fit, vcov = NULL)
        if (termPredictors) fit <- c(fit, termPredictors = NULL)
    }
    else {
        onlyLin <- identical(attr(modelTerms, "prefixLabels"), "")
        if (onlyLin) {
            X <- model.matrix(modelTerms, modelData)
            coefNames <- colnames(X)
        }
        else {
            modelTools <- gnmTools(modelTerms, modelData,
                                   method == "model.matrix", termPredictors)
            coefNames <- names(modelTools$classID)
        }
        if (method == "coefNames") return(coefNames)
        nParam <- length(coefNames)

        if (is.null(constrain))
          constrain <- rep.int(FALSE, nParam)
        else {
            if (identical(constrain, "pick")) {
                if (!require(tcltk))
                    stop("constrain = \"pick\", and tcltk not installed")
                if (!require(relimp))
                    stop("the relimp package from CRAN needs to be installed")
                call$constrain <-
                    pickFrom(coefNames[seq(coefNames) > nElim],
                             setlabels = "Coefficients to constrain",
                             title = "Constrain one or more gnm coefficients",
                             items.label = "Model coefficients:",
                             edit.setlabels = FALSE)
                call$constrain <- unname(unlist(call$constrain))
                if(!length(nchar(call$constrain))) {
                    warning("no parameters were specified to constrain")
                    call$constrain <- NULL
                }
                constrain <- is.element(coefNames, call$constrain)
            }
            else if (is.character(constrain)) {
                constrain <- is.element(coefNames, constrain)
            }
            else if (is.numeric(constrain)) {
                if (!missing(eliminate) & any(constrain < nElim))
                    stop("'constrain' specifies one or more parameters",
                         "in 'eliminate' term(s)")
                constrain <- is.element(seq(nParam), constrain)
            }
            if (length(constrain) != nParam)
                stop("length of 'constrain' not equal to number of parameters")
        }

        if (is.null(start))
            start <- rep.int(NA, nParam)
        else if (length(start) != nParam) {
            if (!missing(eliminate) & length(start) == (nParam - nElim))
                start <- c(rep(NA, nElim), start)
            else
                stop("length(start) must either equal the no. of parameters\n",
                     "or the no. of non-eliminated parameters.")
        }

        if (onlyLin)
            X[, constrain] <- 0
        else {
            theta <- seq(start)
            theta[!is.na(start)] <- start[!is.na(start)]
            theta[constrain] <- 0
            factorList <- modelTools$factorList(theta)
            X <- modelTools$localDesignFunction(theta, factorList)
            attr(X, "assign") <- modelTools$termAssign
        }
        if (method == "model.matrix") return(X)

        if (onlyLin) {
            if (any(is.na(start))) start <- NULL
            fit <- glm.fit(X, y, family = family, weights = weights,
                           offset = offset, start = start,
                           control = do.call("glm.control",
                           unname(gnmControl(...)[-2])),
                           intercept = attr(attr(modelTerms, "terms"),
                           "intercept"))
            if (x) fit$x <- X
            if (vcov) fit$vcov <- stats:::vcov.glm(fit)
            fit <- fit[-c(4,5,7,12,17,20)]
            names(fit)[6] <- "predictors"
        }
        else
            fit <- gnmFit(modelTools, y, constrain, nElim, family, weights,
                          offset, nObs = nObs, start = start,
                          control = gnmControl(...), verbose, x, vcov,
                          termPredictors)
    }
    if (is.null(fit)) {
        warning("Algorithm failed - no model could be estimated", call. = FALSE)
        return()
    }
    fit <- c(list(call = call, formula = formula,
                  terms = attr(modelTerms, "terms"), constrain = constrain,
                  eliminate = nElim,  na.action = attr(modelData, "na.action"),
                  xlevels = .getXlevels(attr(modelData, "terms"), modelData),
                  offset = offset, control = control), fit)

    asY <- c("predictors", "fitted.values", "residuals", "prior.weights",
             "weights", "y", "offset")
    if (inherits(data, "table") & !is.empty.model(modelTerms)) {
        fit[asY] <- lapply(fit[asY], replace,
                           list = as.numeric(names(fit$y)), x = data)
        if (!is.null(fit$na.action)) fit$na.action <- NULL
    }
    else
        fit[asY] <- lapply(fit[asY], structure, names = names(y))
    
    if (model) {
        attr(modelData, "terms") <- attr(modelTerms, "terms")
        fit$model <- modelData
    }
    class(fit) <- c("gnm", "glm", "lm")
    fit
}

