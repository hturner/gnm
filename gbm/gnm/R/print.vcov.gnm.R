print.vcov.gnm <- function(x, ...) {
    if (attr(x, "eliminate")) {
        if (attr(x, "eliminate") < nrow(x)) { 
            keep <- (attr(x, "eliminate") + 1):nrow(x)
            print.default(x[keep, keep])
        }
        else
            cat("No non-eliminated coefficients\n")
    }
    else
        print.default(x)
}
