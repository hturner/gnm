unlistOneLevel <- function(theList){
     result <- vector(length = sum(sapply(theList,
                      function(x) if(is.list(x)) length(x) else 1)),
                      mode = "list")
     count <- 0
     for (i in seq(theList)){
         theItem <- theList[[i]]
         if (is.list(theItem)){
             for (j in seq(theItem)){
                 count <- count + 1
                 result[[count]] <- theItem[[j]]
             }
         }
         else {
             count <- count + 1
             result[[count]] <- theItem
         }
     }
     return(result[1:count])
}
