cautres <- read.table("cautres.txt")
cautres <-
    as.table(tapply(cautres$count,
                    cautres[c("vote", "class", "religion", "election")], sum))
