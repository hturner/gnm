print.vcov.gnm <- function(x, ...) {
  print.default(x[!attr(x, "auxiliary"),
                       !attr(x, "auxiliary") ])
}
