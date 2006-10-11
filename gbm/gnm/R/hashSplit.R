hashSplit <- function(string){
    ## An adaptation of some Python code by 'tim'
    ## http://forum.textdrive.com/viewtopic.php?id=3095
    if (!length(string) || !nchar(string))
        return(string)
    s <- strsplit(string, "")[[1]]
    a <- 0
    ans <- vector("list", length(s))
    iq  <- FALSE
    for (z in seq(s)) {
        if (s[z] ==  "#" & !iq) {
            ans[z] <- paste(s[a:(z - 1)], collapse = "")
            a <- z + 1
        }
        else if (s[z] == "\""){
            iq <- !iq
        }
    }
    ans[z] <- paste(s[a:z], collapse = "")
    unlist(ans)
}
