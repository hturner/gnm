print.vcov.gnm <- function(x, ...) {
    if (!is.null(attr(x, "ofInterest")))
        print.default(x[attr(x, "ofInterest"), attr(x, "ofInterest")])
    else
        print.default(x)
}
