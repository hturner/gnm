print.summary.gnm <- function (x, digits = max(3, getOption("digits") - 3),
                               signif.stars = getOption("show.signif.stars"),
                               symbolic.cor = x$symbolic.cor, ...)
{
    cat("\nCall:\n", deparse(x$call), "\n", sep = "", fill = TRUE)

    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q",
            "Max")
    }
    print.default(x$deviance.resid, digits = digits, na = "", print.gap = 2)

    tidy.zeros <- function(vec)
        ifelse(abs(vec) < 100 * .Machine$double.eps, 0, vec)
    coefs <- tidy.zeros(coef(x))
    if (!is.null(ofInterest(x)))
        coefs <- coefs[ofInterest(x), , drop = FALSE]
    non.elim <- nrow(coefs)
    elim <- length(x$elim.coefs)

    if (non.elim | elim) {
        cat("\nCoefficients", " of interest"[!is.null(ofInterest(x))], ":\n",
            sep = "")
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
            signif.legend = !elim, na.print = "NA", ...)
        if (elim){
            cat("\nEliminated coefficients:\n", sep = "")
            printCoefmat(x$elim.coefs, digits = digits,
                         signif.stars = signif.stars, na.print = "NA", ...)
        }
        coefs <- c(coefs[,2], x$elim.coefs[,2])
        if (any(!is.na(coefs)))
            cat("\n(Dispersion parameter for ", x$family$family,
                " family taken to be ", format(x$dispersion), ")\n", sep = "")
        if (any(is.na(coefs)))
            cat("\nStd. Error is NA where coefficient has been constrained or",
                "is unidentified\n")
    }
    else cat("\nNo coefficients", " of interest"[!is.null(ofInterest(x))],
             ". \n\n", sep = "")

    cat("\nResidual deviance: ", format(x$deviance,
                                        digits = max(5, digits + 1)),
        " on ", format(x$df.residual, digits = max(5, digits + 1)),
        " degrees of freedom\n",
        "AIC: ", format(x$aic, digits = max(4, digits + 1)), "\n\n",
        "Number of iterations: ", x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        if (attr(x$cov.scaled, "eliminate")) {
            eliminate <- seq(attr(x$cov.scaled, "eliminate"))
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
