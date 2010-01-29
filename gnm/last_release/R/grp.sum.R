grp.sum <- function(x, grp.end){
    x <- cumsum(x)[grp.end]
    x - c(0, x[-length(x)])
}
