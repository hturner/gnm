getData <- function() {
    nFrame <- match(TRUE, sapply(sys.calls(), function(x) {
        identical(x[[1]], as.name("gnmTerms"))}))  
    get("data", sys.frame(nFrame))
}
    
