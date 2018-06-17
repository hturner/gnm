#  Copyright (C) 2005, 2006, 2010 David Firth and Heather Turner
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

## now only computes se for non-eliminated parameters
se <- function(object, ...) {
    UseMethod("se", object)
}

se.default <- function(object, ...){
    stop("No se method defined for this class of object")
}

se.gnm <- function(object, estimate = NULL, checkEstimability = TRUE,
                   Vcov = NULL, dispersion = NULL, ...){
    if (!is.null(Vcov) && !is.null(dispersion)){
        Vcov <- Vcov * dispersion
    } else {
        Vcov <- vcov(object, dispersion = dispersion, use.eliminate = FALSE)
    }
    if (!length(Vcov)) return("Model has no non-eliminated parameters")
    coefs <- coef(object)
    coefNames <- names(coefs)
    eliminate <- object$eliminate
    nelim <- nlevels(eliminate)
    l <- length(coefs)
    if (identical(estimate, "[?]"))
        estimate <- pickCoef(object,
                             title = paste("Estimate standard errors",
                             "for one or more gnm coefficients"))
    if (is.null(estimate)){
        if (!is.null(object$ofInterest)) estimate <- ofInterest(object)
        else estimate <- seq(object$coefficients)
    }
    if (is.character(estimate))
        estimate <- match(estimate, coefNames, 0)
    if (is.vector(estimate) && all(estimate %in% seq(coefs))) {
        if (!length(estimate))
            stop("no non-eliminated coefficients specified by 'estimate'",
                 "argument")
        comb <- naToZero(coefs[estimate])
        var <- Vcov[estimate, estimate]
        coefMatrix <- matrix(0, l, length(comb))
        coefMatrix[cbind(estimate, seq(length(comb)))] <- 1
        colnames(coefMatrix) <- names(comb)
    }
    else {
        coefMatrix <- as.matrix(estimate)
        if (!is.numeric(coefMatrix))
            stop("'estimate' should specify parameters using ",
                 "\"pick\" or a vector of \n names/indices; ",
                 "or specify linear combinations using ",
                 "a numeric vector/matrix.")
        if (nrow(coefMatrix) != l)
            stop("NROW(estimate) should equal ",
                 "length(coef(model)) - nlevels(model$eliminate)")
        comb <- drop(crossprod(coefMatrix, naToZero(coefs)))
        var <- crossprod(coefMatrix, crossprod(Vcov, coefMatrix))
    }
    estimable <- rep(TRUE, ncol(coefMatrix))
    if (checkEstimability) {
        estimable <- checkEstimable(object, coefMatrix, ...)
        if (any(!na.omit(estimable)))
            message("Std. Error is NA where estimate is fixed or ", 
                    "unidentified")
    }
    if (is.matrix(var))
        sterr <- sqrt(diag(var))
    else
        sterr <- sqrt(var)
    is.na(sterr[estimable %in% c(FALSE, NA)]) <- TRUE
    result <- data.frame(comb, sterr)
    rowNames <- colnames(coefMatrix)
    if (is.null(rowNames))
        rowNames <- paste("Combination", ncol(coefMatrix))
    dimnames(result) <- list(rowNames, c("Estimate", "Std. Error"))
    result
}
