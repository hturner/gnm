print.coef.gnm <- function(x, ...) {
    cat("Coefficients", " of interest"[!is.null(attr(x, "ofInterest"))], ":\n",
            sep = "")
    if (!is.null(attr(x, "ofInterest")))
        print.default(format(x[attr(x, "ofInterest")]), quote = FALSE)
    else
        print.default(format(x), quote = FALSE)
}
