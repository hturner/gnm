print.coef.gnm <- function(x, ...) {
    print.default(x[!attr(x, "auxiliary")])
}
