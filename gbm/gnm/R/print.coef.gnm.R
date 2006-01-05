print.coef.gnm <- function(x, ...) {
    if (attr(x, "eliminate"))
        print.default(x[(attr(x, "eliminate") + 1):length(x)])
    else
        print.default(x)
}
