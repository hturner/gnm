print.vcov.gnm <- function(object) {
  print.default(object[!attr(object, "auxiliary"),
                       !attr(object, "auxiliary") ])
}
