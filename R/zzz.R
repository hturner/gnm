.onUnload <- function(libpath) {
    library.dynam.unload("gnm", libpath)
}

messageVector <- function(x){
    message(paste(strwrap(paste(x, collapse = ", ")), 
                  collapse = "\n"))
}
