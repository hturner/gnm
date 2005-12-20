print.coef.gnm <- function(x, ...) {
    if (attr(x, "eliminate"))
        print.default(x[-seq(attr(x, "eliminate"))])
    else
        print.default(x)
}
