mentalHealth <- read.table("mentalHealth.txt")
mentalHealth$SES <- as.ordered(mentalHealth$SES)
mentalHealth$MHS <- ordered(mentalHealth$MHS,
                            levels=c("well", "mild", "moderate", "impaired"))

                               
