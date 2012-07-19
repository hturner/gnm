# this should always be a summary based on single grouping factor
summary.meanResiduals <- function (object, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nModel call:\n", deparse(attr(object, "call"),
                                   width.cutoff = options()$width),
        "\n", sep = "", fill = TRUE)

    cat("Mean residuals by ", attr(object, "by"),  ":\n\n", sep = "")
    q <- quantile(object, na.rm = TRUE)
    names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
    print.default(q, digits = digits, na.print = "", print.gap = 2)

    if (attr(object, "standardized")) {
        cat("\nTest of Normality:\n")
        df <- attr(object, "df")
        if (df > 0) {
            chi.sq <- sum(as.vector(object)^2)
            p.value <- pchisq(chi.sq, df, lower.tail = FALSE)
            test <- c(chi.sq, df, p.value)
            cat("\nChi^2 =", format(chi.sq, digits = digits), "on",
                df, "df, p-value =", format(p.value, digits = digits), "\n")
        }
        else cat("\n(zero degrees of freedom)\n")
    }
    else cat("\nResiduals are not standardized\n")
}
