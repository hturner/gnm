print.vcov.gnm <- function(x, ...) {
    if (attr(x, "eliminate")) {
        eliminate <- seq(attr(x, "eliminate"))
        print.default(x[-eliminate, -eliminate])
    }
    else
        print.default(x)
}
