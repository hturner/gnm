quick.rowsum <- function (x, group, ugroup,...)
{
    .Call("Rrowsum_matrix", x, NCOL(x), group, ugroup,
          FALSE, PACKAGE = "base")
}
