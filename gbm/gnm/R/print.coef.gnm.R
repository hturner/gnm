print.coef.gnm <- function(x, ...) {
    cat("Coefficients", " of interest"[!is.null(ofInterest(x))], ":\n",
            sep = "")
    if (!is.null(ofInterest(x)))
        print.default(format(x[ofInterest(x)]), quote = FALSE)
    else
        print.default(format(x), quote = FALSE)
}
