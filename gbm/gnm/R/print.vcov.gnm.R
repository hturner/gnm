print.vcov.gnm <- function(x, ...) {
    if (attr(x, "eliminate")) {
        keep <- (attr(x, "eliminate") + 1):nrow(x)
        print.default(x[-eliminate, -eliminate])
    }
    else
        print.default(x)
}
