meanResiduals <- function(object, by = NULL, standardized = TRUE,
                          as.table = TRUE, ...){
    if (is.null(by))
        stop("`by' must be specified in order to compute grouped residuals")
    if (inherits(by, "formula")){
        ## check single factor only
        if (ncol(attr(terms(by), "factors")) != 1)
            stop("`by' should only specify a single term")
        ## find factors as in mosaic.glm
        by <- do.call("model.frame",
                      list(formula = by, data = object$data,
                           subset = object$call$subset,
                           na.action = na.pass, drop.unused.levels = TRUE))
        ## following loop needed due to bug in model.frame.default (fixed for R 2.12)
        for(nm in names(by)) {
            f <- by[[nm]]
            if(is.factor(f) && length(unique(f[!is.na(f)])) < length(levels(f)))
                by[[nm]] <- by[[nm]][, drop = TRUE]
        }
        if (!is.null(object$na.action))
            by <- by[-object$na.action,]
    }
    if (!all(sapply(by, is.factor)))
        warning("Coercing variables specified by `by' to factors")
    fac <- factor(interaction(by)) # drop unused levels
    if (length(fac) != length(object$y))
        stop("Grouping factor of length", length(fac),
             "but model frame of length", length(object$y))

    r <- object$residuals
    ## recompute weights for better accuracy
    w  <- as.numeric(object$prior.weights *
                     object$family$mu.eta(predict(object, type = "link"))^2/
                     object$family$variance(object$fitted))
    agg.wts <- tapply(w, by, sum) #unlike rowsum, keeps all levels of interaction
    res <- tapply(r * w, by, sum)/agg.wts
    if (standardized) res <- res * sqrt(agg.wts)
    ## now compute degrees of freedom
    Xreduced <- rowsum(model.matrix(object), fac, na.rm = TRUE)
    ## suppressWarnings in rankMatrix re coercion to dense matrix
    if (as.table){
        res <- structure(as.table(res), call = object$call,
                         by = paste(names(by), collapse = ":"),
                         df = nlevels(fac) -
                         suppressWarnings(rankMatrix(Xreduced)),
                         standardized = standardized,
                         weights = as.table(agg.wts))
        class(res) <- c("meanResiduals", "table")
    }
    else {
        res <- structure(c(res), call = object$call,
                         by = paste(names(by), collapse = ":"),
                         df = nlevels(fac) -
                         suppressWarnings(rankMatrix(Xreduced)),
                         standardized = standardized,
                         weights = c(agg.wts))
        class(res) <- c("meanResiduals", "numeric")
    }
    return(res)
}



