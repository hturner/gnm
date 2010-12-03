print.meanResiduals <- function (object, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nModel call:\n", deparse(object$call), "\n", sep = "", fill = TRUE)

    cat("Mean residuals by", deparse(object$by),  ":\n\n")
    if (length(object$residuals > 20)) {
        object$residuals <- quantile(object$residuals, na.rm = TRUE)
        names(object$residuals) <- c("Min", "1Q", "Median", "3Q",
                                     "Max")
    }
    print.default(object$residuals, digits = digits, na = "", print.gap = 2)

    cat("\nTest of Normality:\n")
    if (object$df > 0) {
        chi.sq <- sum(as.vector(object$residuals)^2)
        p.value <- pchisq(chi.sq, object$df, lower.tail = FALSE)
        test <- c(chi.sq, object$df, p.value)
        cat("\nChi^2 =", format(chi.sq, digits = digits), "on",
            object$df, "df, p-value =", format(p.value, digits = digits), "\n")
    }
    else cat("\n(zero degrees of freedom)\n")
}
