DrefVariables <- function(..., formula = ~ 1) {
    c(match.call(expand.dots = FALSE)[[2]],
      as.list(attr(terms(formula), "variables"))[-1])
}
    
