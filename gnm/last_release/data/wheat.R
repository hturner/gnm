wheat <- read.table("wheat.txt")
wheat$year <- ordered(wheat$year)
wheat$N <- ordered(wheat$N, levels = c("0", "n", "N"))

