backpain <- read.table("backpain.txt")
backpain$pain <- ordered(backpain$pain,
                         levels = c("worse", "same", "slight.improvement",
                         "moderate.improvement", "marked.improvement",
                         "complete.relief"))
