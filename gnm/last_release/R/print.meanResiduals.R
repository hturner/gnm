# this should just print residuals as default (could be as table)
print.meanResiduals <- function (object, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nModel call:\n", deparse(object$call, width = options()$width),
        "\n", sep = "", fill = TRUE)

    cat("Mean residuals by ", object$by,  ":\n\n", sep = "")
    print(object$residuals, digits = digits, print.gap = 2)
}
