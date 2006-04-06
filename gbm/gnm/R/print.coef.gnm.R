print.coef.gnm <- function(x, ...) {
    if (attr(x, "eliminate")) {
        if (attr(x, "eliminate") < length(x))
            print.default(x[(attr(x, "eliminate") + 1):length(x)])
        else
            cat("No non-eliminated coefficients\n")
    }
    else {
        print.default(format(x), quote = FALSE)
    }
}
