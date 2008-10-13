Const <- function(const) {
    badCall <- !"nonlinTerms" %in% lapply(sys.calls(), "[[", 1)
    if (any(badCall))
        stop("Const terms are only valid in the predictors of \"nonlin\" ",
             "functions.")
    call("rep", substitute(const), quote(nObs))
}

