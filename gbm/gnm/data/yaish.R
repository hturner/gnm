yaish <- read.table("yaish.txt")
yaish <- as.table(tapply(yaish$count,
                         yaish[c("educ", "orig", "dest")], sum))
