print.coef.gnm <- function(x, ...) {
    if (!is.null(attr(x, "eliminate")))
        print.default(x[-seq(attr(x, "eliminate"))])
    else
        print.default(x)
}
