occupational.status <- read.table("occupational.status.txt")
occupational.status <-
    as.table(tapply(occupational.status$counts,
                    occupational.status[c("origin", "destination")], sum))

