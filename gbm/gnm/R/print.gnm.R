print.gnm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n", deparse(x$call), "\n", sep = "", fill = TRUE)    
    if (length(ofInterest(x)) || length(coef(x))) {
        cat("Coefficients", " of interest"[!is.null(ofInterest(x))], ":\n",
            sep = "")
        if (!is.null(ofInterest(x)))
            print.default(format(coef(x)[ofInterest(x)], digits = digits),
                          print.gap = 2, quote = FALSE)
        else
            print.default(format(coef(x), digits = digits), print.gap = 2,
                          quote = FALSE)
    }
    else cat("No coefficients. \n\n", sep = "")
    cat("\nDeviance:           ", format(x$deviance, digits),
        "\nPearson chi-squared:",
        format(sum(na.omit(c(residuals(x, type = "pearson")))^2), digits),
        "\nResidual df:        ", x$df.residual, "\n")
    invisible(x)
}
