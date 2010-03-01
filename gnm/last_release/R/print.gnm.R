print.gnm <- function (x, digits = max(3, getOption("digits") - 3),
                       show.elim = FALSE, ...) {
    cat("\nCall:\n", deparse(x$call), "\n", sep = "", fill = TRUE)
    non.elim <- length(coef(x)) && length(ofInterest(x))
    elim <- show.elim && length(x$elim.coefs)
    if (non.elim) {
        cat("Coefficients", " of interest"[!is.null(ofInterest(x))], ":\n",
            sep = "")
        if (!is.null(ofInterest(x)))
            print.default(format(coef(x)[ofInterest(x)], digits = digits),
                          print.gap = 2, quote = FALSE)
        else
            print.default(format(coef(x), digits = digits), print.gap = 2,
                          quote = FALSE)

    }
    if (elim){
        cat("\nEliminated coefficients:\n")
        print.default(format(x$elim.coefs, digits = digits), print.gap = 2,
                      quote = FALSE)
    }
    if (!non.elim && !elim)
        cat("No coefficients", " of interest"[!is.null(ofInterest(x))],
             ". \n\n", sep = "")
    cat("\nDeviance:           ", format(x$deviance, digits),
        "\nPearson chi-squared:",
        format(sum(na.omit(c(residuals(x, type = "pearson")))^2), digits),
        "\nResidual df:        ", x$df.residual, "\n")
    invisible(x)
}
