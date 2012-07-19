# this should just print residuals as default (could be as table)
print.meanResiduals <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nModel call:\n", deparse(attr(x, "call"),
                                   width.cutoff = options()$width),
        "\n", sep = "", fill = TRUE)

    cat("Mean residuals by ", attr(x, "by"),  ":\n\n", sep = "")
    if (!inherits(x, "table")) x <- as.numeric(x)
    NextMethod(object = x, digits = digits, print.gap = 2, ...)
}
