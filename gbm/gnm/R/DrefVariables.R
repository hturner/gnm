DrefVariables <- function(..., formula = ~ 1) {
    c(as.character((match.call(expand.dots = FALSE))[[2]]),
      rownames(attr(terms(formula), "factors")))
}
    
