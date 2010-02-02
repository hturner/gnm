print.coef.gnm <- function(x, ...) {
    if (!is.null(attr(x, "ofInterest"))) {
        if (length(attr(x, "ofInterest"))){
            cat("Coefficient of interest:\n", sep = "")
            print.default(format(x[attr(x, "ofInterest")]), quote = FALSE)
        }
        else
            cat("No coefficients of interest\n")
    }
    else {
        cat("Coefficients:\n")
        print.default(format(x), quote = FALSE)
    }
}
