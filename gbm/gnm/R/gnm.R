gnm <- function(formula, eliminate = NULL, ofInterest = NULL,
                constrain = numeric(0),
                constrainTo = numeric(length(constrain)), family = gaussian,
                data = NULL, subset, weights, na.action,  method = "gnmFit",
                offset, start = NULL, tolerance = 1e-6, iterStart = 2,
                iterMax = 500, trace = FALSE, verbose = TRUE, model = TRUE,
                x = FALSE, termPredictors = FALSE, lsMethod = "qr", ...) {

    call <- match.call()

    modelTerms <- gnmTerms(formula, substitute(eliminate), data)
    modelData <- match.call(expand.dots = FALSE)
    argPos <- match(c("data", "subset", "weights", "na.action", "offset"),
                    names(modelData), 0)
    modelData <- as.call(c(as.name("model.frame"),
                           formula = modelTerms,
                           as.list(modelData)[argPos],
                           drop.unused.levels = TRUE))
    if (inherits(data, "table") && !is.empty.model(modelTerms)) {
        xFactors <- modelData
        xFactors$formula <- Freq ~ .
        xFactors <- eval(xFactors, parent.frame())[, -1]
    }
    modelData <- eval(modelData, parent.frame())

    if (!missing(eliminate)) {
        Elim <- suppressWarnings(eval(substitute(eliminate), modelData))
        if (!is.factor(Elim))
            stop("'eliminate' must be a factor")
        nElim <- nlevels(Elim)
        if (missing(lsMethod)) lsMethod <- "chol"
    }
    else nElim <- 0

    if (method == "model.frame") {
        attr(modelData, "terms") <- attr(modelTerms, "terms")
        return(modelData)
    }
    else if (!method %in% c("gnmFit", "coefNames", "model.matrix") &&
             !is.function(get(method))) {
        warning("function ", method, " can not be found. Using \"gnmFit\".\n",
                call. = FALSE)
        method <- "gnmFit"
    }

    y <- model.response(modelData, "numeric")
    nObs <- NROW(y)

    weights <- as.vector(model.weights(modelData))
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights are not allowed")
    if (is.null(weights))
        weights <- rep.int(1, nObs)
    offset <- as.vector(model.offset(modelData))
    if (is.null(offset))
        offset <- rep.int(0, nObs)

    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }

    if (family$family == "binomial") {
        if (is.factor(y) && NCOL(y) == 1)
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
            return(matrix(, nrow(modelData), 0))
        if (!family$valideta(offset))
            stop("invalid predictor values in empty model")
        mu <- family$linkinv(offset)
        if (!family$validmu(mu))
            stop("invalid fitted values in empty model")
        dmu <- family$mu.eta(offset)
        dev <- sum(family$dev.resids(y, mu, weights))
        modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nObs), mu,
                                                weights, dev))
        fit <- list(coefficients = numeric(0), constrain = numeric(0),
                    constrainTo = numeric(0), eliminate = 0,
                    predictors = offset, fitted.values = mu, deviance = dev,
                    aic = modelAIC, iter = 0, conv = NULL,
                    weights = weights*dmu^2/family$variance(mu),
                    residuals = (y - mu)/dmu, df.residual = nObs, rank = 0,
                    family = family, prior.weights = weights, y = y,
                    converged = NA)
        if (x) fit <- c(fit, x = matrix(, nrow(modelData), 0))
        if (termPredictors) fit <- c(fit, termPredictors =
                                     matrix(, nrow(modelData), 0))
    }
    else {
        onlyLin <- nElim == 0 && all(attr(modelTerms, "prefixLabels") == "")
        if (onlyLin) {
            X <- model.matrix(attr(modelTerms, "terms"), modelData)
            coefNames <- colnames(X)
        }
        else {
            modelTools <- gnmTools(modelTerms, modelData,
                                   method == "model.matrix", termPredictors)
            coefNames <- names(modelTools$classID)
        }
        if (method == "coefNames") return(coefNames)
        nParam <- length(coefNames)

        if (identical(constrain, "[?]"))
            call$constrain <- constrain <- 
                relimp:::pickFrom(coefNames,
                                  subset = (nElim + 1):nParam,
                                  setlabels = "Coefficients to constrain",
                                  title =
                                  "Constrain one or more gnm coefficients",
                                  items.label = "Model coefficients:",
                                  warningTest =
                                  "No parameters were specified to constrain",
                                  return.indices = TRUE)
        if (is.character(constrain)) {
            if (length(constrain) == 1)
                constrain <- grep(constrain, coefNames)
            else
                constrain <- match(constrain, coefNames)
        }
        ## dropped logical option
        if (any(constrain < nElim))
            stop("'constrain' specifies one or more parameters",
                 "in the 'eliminate' term.")
        if (!all(constrain %in% seq(coefNames)))
            stop(" 'constrain' specifies non-existant parameters.")
        
        if (is.null(start))
            start <- rep.int(NA, nParam)
        else if (length(start) != nParam) {
            if (!missing(eliminate) && length(start) == (nParam - nElim))
                start <- c(rep(NA, nElim), start)
            else
                stop("length(start) must either equal the no. of parameters\n",
                     "or the no. of non-eliminated parameters.")
        }
        
        if (onlyLin) {
            offset <- offset + X[, constrain, drop = FALSE] %*% constrainTo
            X[, constrain] <- 0
        }
        else {
            theta <- seq(start)
            theta[!is.na(start)] <- start[!is.na(start)]
            theta[constrain] <- constrainTo
            factorList <- modelTools$factorList(theta)
            X <- modelTools$localDesignFunction(theta, factorList)
            attr(X, "assign") <- modelTools$termAssign
        }
        if (method == "model.matrix")
            return(X)

        if (!is.numeric(tolerance) || tolerance <= 0)
            stop("value of 'tolerance' must be > 0")
        if (!is.numeric(iterMax) || iterMax <= 0)
            stop("maximum number of iterations must be > 0")

        if (onlyLin) {
            if (any(is.na(start))) start <- NULL
            if (verbose) cat("Linear predictor - using glm.fit\n")
            fit <- glm.fit(X, y, family = family, weights = weights,
                           offset = offset, start = start,
                           control = glm.control(tolerance, iterMax, trace),
                           intercept = attr(attr(modelTerms, "terms"),
                           "intercept"))
            if (sum(is.na(coef(fit))) > length(constrain)) {
                extra <- setdiff(which(is.na(coef(fit))), constrain)
                ind <- order(c(constrain, extra))
                constrain <- c(constrain, extra)[ind]
                constrainTo <- c(constrainTo, numeric(length(extra)))[ind]
            }
            fit$constrain <- constrain
            fit$constrainTo <- constrainTo
            if (x) fit$x <- X
            fit <- fit[-c(4,5,7,12,17,20)]
            names(fit)[6] <- "predictors"
        }
        else if (method != "gnmFit")
            fit <- do.call(method, list(modelTools, y, constrain, constrainTo,
                                        nElim, family, weights, offset, nObs,
                                        start, tolerance, iterStart, iterMax,
                                        trace, verbose, x, termPredictors, ...))
        else
            fit <- gnmFit(modelTools, y, constrain, constrainTo, nElim, family,
                          weights, offset, nObs, start, tolerance, iterStart,
                          iterMax, trace, verbose, x, termPredictors,
                          lsMethod = lsMethod)
    }
    if (is.null(fit)) {
        warning("Algorithm failed - no model could be estimated", call. = FALSE)
        return()
    }
    
    if (is.null(ofInterest) && !missing(eliminate))
        ofInterest <- (nElim + 1):length(coefNames)
    if (identical(ofInterest, "[?]"))
        call$ofInterest <- ofInterest <- 
            pickCoef(coefNames, 
                     warningText = paste("No subset of coefficients selected",
                     "- assuming all are of interest."))
    if (is.character(ofInterest)) {
        if (length(ofInterest) == 1)
            ofInterest <- grep(ofInterest, coefNames)
        else
            ofInterest <- match(ofInterest, coefNames)
    }
    if (!is.null(ofInterest)) {
        if (!any(ofInterest %in% seq(coefNames))) 
            stop("'ofInterest' does not specify a subset of the ",
                 "coefficients.")
        names(ofInterest) <- coefNames[ofInterest]
    }

    fit <- c(list(call = call, formula = formula,
                  terms = attr(modelTerms, "terms"), eliminate = nElim,
                  ofInterest = ofInterest,
                  na.action = attr(modelData, "na.action"),
                  xlevels = .getXlevels(attr(modelData, "terms"), modelData),
                  offset = offset, tolerance = tolerance, iterStart = iterStart,
                  iterMax = iterMax), fit)

    asY <- c("predictors", "fitted.values", "residuals", "prior.weights",
             "weights", "y", "offset")
    if (exists("xFactors", inherits = FALSE)) {
        fit[asY] <- lapply(fit[asY], tapply, xFactors, sum)
        if (!is.null(fit$na.action)) fit$na.action <- NULL
    }
    else
        fit[asY] <- lapply(fit[asY], structure, names = names(y))

    if (model) {
        attr(modelData, "terms") <- attr(modelTerms, "terms")
        fit$model <- modelData
    }
    class(fit) <- c("gnm", "glm", "lm")
    attr(fit, ".Environment") <- environment(gnm)
    fit
}

