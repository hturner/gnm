print.summary.gnm <- function (x, digits = max(3, getOption("digits") - 3),
                                symbolic.cor = x$symbolic.cor, ...) 
{
    if (!is.null(x$original.call))
      cat("\nCall:\n", deparse(x$original.call), "\n\n", sep = "",
          fill = TRUE)
    else
      cat("\nCall:\n", deparse(x$call), "\n\n", sep = "", fill = TRUE)
    
    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", 
            "Max")
    }
    print.default(x$deviance.resid, digits = digits, na = "", print.gap = 2)
    
    if (length(coef(x))) {
        cat("\nCoefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2,
                      quote = FALSE)
    }
    else cat("\nNo coefficients\n\n")
    
    cat("\nResidual deviance: ", format(x$deviance,
                                        digits = max(5, digits + 1)),
        " on ", format(x$df.residual, digits = max(5, digits + 1)),
        " degrees of freedom\n",
        "AIC: ", format(x$aic, digits = max(4, digits + 1)), "\n\n",
        "Number of iterations: ", x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
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
