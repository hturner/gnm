print.profile.gnm <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    #if (attr(x, "eliminate"))
       # coefs <- coefs[-seq(attr(x$cov.scaled, "eliminate")), ]
    
    if (length(x)) {
        if (any(sapply(x, is.null)))
            cat("\nProfile is NULL where coefficient has been constrained or",
                "is unidentified\n\n")
        print.default(x)
    }
    else cat("\nNo ", "non-eliminated "[attr(x$cov.scaled, "eliminate") > 0],
             "coefficients\n\n", sep = "")
    invisible(x)
}
