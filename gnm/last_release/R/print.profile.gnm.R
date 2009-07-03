print.profile.gnm <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    #if (attr(x, "eliminate"))
       # coefs <- coefs[-seq(attr(x$cov.scaled, "eliminate")), ]

    if (length(x)) {
        if (any(sapply(x, function(x) isTRUE(is.na(x)))))
            cat("\nProfile is NA where coefficient has been constrained or",
                "is unidentified\n\n")
        print.default(x)
    }
    else cat("\nNo coefficients profiled.\n\n", sep = "")
    invisible(x)
}
