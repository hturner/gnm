getModelFrame <- function() {
    nFrame <- match(TRUE, sapply(sys.calls(),
                       function(x)identical(x[[1]], as.name("gnmTools"))))
    get("gnmData", sys.frame(nFrame))
}
