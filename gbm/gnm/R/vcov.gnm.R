vcov.gnm <-  function(object, ...){
    structure(summary(object, corr = FALSE, ...)$cov.scaled,
              class = "vcov.gnm", auxiliary = object$auxiliary)
}
