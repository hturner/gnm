print.gnm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "", fill = TRUE)    
    if (length(coef(x)) - x$eliminate) {
        cat("Coefficients:\n")
        if (attr(coef(x), "eliminate"))
            print.default(format(coef(x)[-seq(attr(coef(x), "eliminate"))],
                                 digits = digits), print.gap = 2,
                          quote = FALSE)
        else
            print.default(format(coef(x), digits = digits), print.gap = 2,
                          quote = FALSE)
    }
    else cat("No ", "non-eliminated "[x$eliminate > 0], "coefficients\n\n",
             sep = "")
    cat("\nDeviance:           ", format(x$deviance, digits),
        "\nPearson chi-squared:",
        format(sum(na.omit(residuals(x,type = "pearson"))^2), digits),
        "\nResidual df:        ", x$df.residual, "\n")
    invisible(x)
}
