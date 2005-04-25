backPain <- read.table("backpain.txt")
backPain$pain <- ordered(backpain$pain,
                         levels = c("worse", "same", "slight.improvement",
                         "moderate.improvement", "marked.improvement",
                         "complete.relief"))
