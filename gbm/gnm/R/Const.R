Const <- function(const) {
    checkCall()
    call("rep", substitute(const), quote(nObs))
}

