local({
       options(device="x11")
       options(repos="http://cran.uk.r-project.org")
       options(pkgType="source")
       old <- getOption("defaultPackages")
       options(defaultPackages = c(old, "tcltk"))
       options(browser = "/usr/bin/firefox")
       options(SweaveHooks = list(eval = function() {
           options(show.signif.stars = FALSE, width = 90)
       }
               ))
   })

.First <- function() cat("\n   Welcome to David's R !\n\n")
