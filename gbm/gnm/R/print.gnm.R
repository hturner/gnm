print.gnm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "", fill = TRUE)
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(x$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nDeviance:           ", x$deviance,
        "\nPearson chi-squared:", sum(na.omit(residuals(x,
                                                        type = "pearson"))^2),
        "\nResidual df:        ", x$df.residual, "\n")
    invisible(x)
}
