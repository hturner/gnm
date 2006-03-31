confint.profile.gnm <- function (object, parm = names(object),
                                 level = 0.95, ...)  {
    of <- attr(object, "original.fit")
    pnames <- names(coef(of))
    if (is.character(parm)) 
        parm <- match(parm, pnames, nomatch = 0)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(round(100 * a, 1), "%")
    ci <- array(NA, dim = c(length(parm), 2), dimnames = list(pnames[parm], 
        pct))
    cutoff <- qnorm(a)
    std.err <- attr(object, "summary")$coefficients[, "Std. Error"]
    parm <- parm[!is.na(std.err)[parm]]
    for (pm in parm) {
        pro <- object[[pnames[pm]]]   
        if (length(pnames) > 1) 
            sp <- spline(x = pro[, "par.vals"][, pm], y = pro[, 
                1])
        else sp <- spline(x = pro[, "par.vals"], y = pro[, 1])
        est <- approx(sp$y, sp$x, xout = cutoff)$y
        ci[pnames[pm], ] <- ifelse(is.na(est) & attr(pro, "asymptote"),
                                   c(-Inf, Inf), est)
    }
    drop(ci)
}
