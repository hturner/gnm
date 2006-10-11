checkCall <- function(){
    badCall <- lapply(sys.calls(), "[[", 1) %in% c("model.frame.default",
                                                   "model.matrix.default")
    if (any(badCall))
        stop(paste(sys.call(-1)[[1]], "terms are only valid in gnm models."))
}
