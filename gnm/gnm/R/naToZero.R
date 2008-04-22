naToZero <- function(vec){
   vec[is.na(vec)] <- 0
   return(vec)
}
