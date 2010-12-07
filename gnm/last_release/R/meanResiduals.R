##' Computes the mean working residuals from a model fitted using
##' Iterative Weighted Least Squares for each level of a factor or
##' interaction of factors.
##'
##' <details>
##'
##' @title Average Residuals within Factor Levels
##' @param object model object for which \code{object$residuals} gives the
##' working residuals and \code{object$weights} gives the working weights.
##' @param by either a formula specifying one or more factors (recommended)
##' or a list of factors (the elements of which must correspond exactly to
##' observations in the model frame). When more than one factor is specified,
##' their interaction is used to specify the grouping factor.
##' @param standardized logical: if \code{TRUE}, the mean residuals are
##' standardized to be approximately standard normal.
##' @param as.table logical: if \code{TRUE} the result is returned as a table
##' cross-classified by the factors passed to \code{by}.
##' @param ... currently ignored
##' @return An object of class \code{"meanResiduals"}, for which \code{print}
##' and \code{summary} methods are provided. A \code{"meanResiduals"}
##' object is a list containing the following elements:
##' \item{ call }{ the call used to create the model object from which the
##' mean residuals are derived. }
##' \item{ by }{ a label for the grouping factor. }
##' \item{ residuals }{ the mean residuals. }
##' \item{ df }{ the degrees of freedom associated with the mean
##' residuals. }
##' \item{ standardized }{ the \code{standardized} argument. }
##' @author Heather Turner
##' @example examples/meanResiduals.R
meanResiduals <- function(object, by = NULL, standardized = TRUE, as.table = TRUE, ...){
    if (is.null(by))
        stop("`by' must be specified in order to compute grouped residuals")
    if (inherits(by, "formula")){
        ## find factors as in mosaic.glm
        by <- do.call("model.frame",
                      list(formula = by, data = object$data, subset = object$call$subset,
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
    fac <- factor(interaction(by)) # drop unused levels
    if (length(fac) != length(object$y))
        stop("Grouping factor of length", length(fac),
             "but model frame of length", length(object$y))

    r <- object$residuals
    ## recompute weights for better accuracy
    w  <- as.numeric(object$prior.weights * object$family$mu.eta(predict(object, type = "link"))^2/
                     object$family$variance(object$fitted))
    agg.wts <- tapply(w, by, sum) #unlike rowsum, keeps all levels of interaction
    res <- tapply(r * w, by, sum)/agg.wts
    if (standardized) res <- res * sqrt(agg.wts)
    if (!as.table){
        res <- structure(c(res), weights = c(agg.wts))
    }
    else
        res <- structure(as.table(res), weights = as.table(agg.wts))
    ## now compute degrees of freedom
    Xreduced <- rowsum(model.matrix(object), fac, na.rm = TRUE)
    res <- list(call = object$call, by = paste(names(by), collapse = ":"), residuals = res,
                df = nlevels(fac) - rankMatrix(Xreduced), standardized = standardized)
    class(res) <- "meanResiduals"
    return(res)
}



