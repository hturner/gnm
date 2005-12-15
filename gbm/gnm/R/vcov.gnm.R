vcov.gnm <-  function(object, ...){
    structure(summary(object, corr = FALSE, ...)$cov.scaled,
              eliminate = object$eliminate, class = "vcov.gnm")
}
