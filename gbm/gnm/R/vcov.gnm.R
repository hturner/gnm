vcov.gnm <-  function(object, ...){
    summary(object, corr = FALSE, ...)$cov.scaled
}
