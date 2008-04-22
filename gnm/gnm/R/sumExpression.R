sumExpression <- function(exprList) {
    expr <- exprList[[1]]
    for (i in seq(exprList)[-1]) {
        expr <- call("+", expr, exprList[[i]])
    }
    expr
}
