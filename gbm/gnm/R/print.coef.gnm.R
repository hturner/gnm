print.coef.gnm <- function(x, ...) {
if (!is.null(x$auxiliary)) print.default(x[!attr(x, "auxiliary")])
  else print.default(x)
}
