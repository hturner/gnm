# this should just print residuals as default (could be as table)
print.meanResiduals <- function (x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nModel call:\n", deparse(x$call, width = options()$width),
        "\n", sep = "", fill = TRUE)

    cat("Mean residuals by ", x$by,  ":\n\n", sep = "")
    print(x$residuals, digits = digits, print.gap = 2)
}
