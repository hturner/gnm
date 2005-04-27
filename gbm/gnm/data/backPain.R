backPain <- read.table("backPain.txt")
backPain$pain <- ordered(backPain$pain,
                         levels = c("worse", "same", "slight.improvement",
                         "moderate.improvement", "marked.improvement",
                         "complete.relief"))
