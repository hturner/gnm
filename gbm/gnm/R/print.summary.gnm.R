print.summary.gnm <- function (x, digits = max(3, getOption("digits") - 3),
                               signif.stars = getOption("show.signif.stars"),
                               symbolic.cor = x$symbolic.cor, ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "", fill = TRUE)
    
    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", 
            "Max")
    }
    print.default(x$deviance.resid, digits = digits, na = "", print.gap = 2)

    coefs <- coef(x)
    if (x$eliminate)
        coefs <- coefs[-seq(x$eliminate), ]
    
    if (nrow(coefs)) {
        cat("\nCoefficients:\n")
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
        if (any(!is.na(coefs[,2])))
            cat("\n(Dispersion parameter for ", x$family$family,
                " family taken to be ", format(x$dispersion), ")\n", sep = "")
        if (any(is.na(coefs[,2])))
            cat("\nStd. Error is NA where parameter has been constrained or",
                "is unidentified\n")
    }
    else cat("\nNo ", "non-eliminated "[attr(coefs, "eliminate") > 0],
             "coefficients\n\n", sep = "")
    
    cat("\nResidual deviance: ", format(x$deviance,
                                        digits = max(5, digits + 1)),
        " on ", format(x$df.residual, digits = max(5, digits + 1)),
        " degrees of freedom\n",
        "AIC: ", format(x$aic, digits = max(4, digits + 1)), "\n\n",
        "Number of iterations: ", x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        if (attr(coefs, "eliminate")) {
            eliminate <- seq(attr(coefs, "eliminate"))
            correl <- correl[-eliminate, -eliminate]
        }
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.col = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}
