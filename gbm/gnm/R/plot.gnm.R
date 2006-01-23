plot.gnm <- function (x, which = c(1:3, 5), caption = c("Residuals vs Fitted", 
    "Normal Q-Q", "Scale-Location", "Cook's distance", "Residuals vs Leverage", 
    "Cook's distance vs Leverage"), panel = if (add.smooth) panel.smooth else points, 
    sub.caption = NULL, main = "", ask = prod(par("mfcol")) < 
        length(which) && dev.interactive(), ..., id.n = 3, labels.id = names(residuals(x)), 
    cex.id = 0.75, qqline = TRUE, cook.levels = c(0.5, 1), add.smooth = getOption("add.smooth"), 
    label.pos = c(4, 2)) 
{
    if (!is.numeric(which) || any(which < 1) || any(which > 5)) 
        stop("'which' must be in 1:5")
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    r <- residuals(x)
    yh <- naresid(x$na.action, na.omit(x$predictors))
    w <- weights(x)
    if (!is.null(w)) {
        wind <- w != 0
        r <- r[wind]
        yh <- yh[wind]
        w <- w[wind]
        labels.id <- labels.id[wind]
    }
    n <- length(r)
    if (any(show[2:6])) {
        s <- sqrt(deviance(x)/df.residual(x))
        hii <- hatvalues(x)
        if (any(show[4:6])) {
            cook <- cooks.distance(x)
        }
    }
    if (any(show[c(2:3, 5)])) {
        ylab23 <- "Std. deviance resid."
        r.w <- if (is.null(w)) 
            r
        else sqrt(w) * r
        rs <- r.w/(s * sqrt(1 - hii))
        rs[is.infinite(rs)] <- NaN
    }
    if (any(show[5:6])) {
        hatval <- hatvalues(x)
        r.hat <- range(hatval, na.rm = TRUE)
        isConst.hat <- diff(r.hat) < 1e-10 * mean(hatval)
    }
    if (any(show[c(1, 3)])) 
        l.fit <- "Predicted values"
    if (is.null(id.n)) 
        id.n <- 0
    else {
        id.n <- as.integer(id.n)
        if (id.n < 0 || id.n > n) 
            stop(gettextf("'id.n' must be in {1,..,%d}", n), 
                domain = NA)
    }
    if (id.n > 0) {
        if (is.null(labels.id)) 
            labels.id <- paste(1:n)
        iid <- 1:id.n
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        if (any(show[2:3])) 
            show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
        text.id <- function(x, y, ind, adj.x = TRUE) {
            labpos <- if (adj.x) 
                label.pos[1 + as.numeric(x > mean(range(x)))]
            else 3
            text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
                pos = labpos, offset = 0.25)
        }
    }
    if (is.null(sub.caption)) {
        cal <- x$call
        if (!is.na(m.f <- match("formula", names(cal)))) {
            cal <- cal[c(1, m.f)]
            names(cal)[2] <- ""
        }
        cc <- deparse(cal, 80)
        nc <- nchar(cc[1])
        abbr <- length(cc) > 1 || nc > 75
        sub.caption <- if (abbr) 
            paste(substr(cc[1], 1, min(75, nc)), "...")
        else cc[1]
    }
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        ylim <- range(r, na.rm = TRUE)
        if (id.n > 0) 
            ylim <- extendrange(r = ylim, f = 0.08)
        plot(yh, r, xlab = l.fit, ylab = "Residuals", main = main, 
            ylim = ylim, type = "n", ...)
        panel(yh, r, ...)
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(caption[1], 3, 0.25)
        if (id.n > 0) {
            y.id <- r[show.r]
            y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
            text.id(yh[show.r], y.id, show.r)
        }
        abline(h = 0, lty = 3, col = "gray")
    }
    if (show[2]) {
        ylim <- range(rs, na.rm = TRUE)
        ylim[2] <- ylim[2] + diff(ylim) * 0.075
        qq <- qqnorm(rs, main = main, ylab = ylab23, ylim = ylim, 
            ...)
        if (qqline) 
            qqline(rs, lty = 3, col = "gray50")
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(caption[2], 3, 0.25)
        if (id.n > 0) 
            text.id(qq$x[show.rs], qq$y[show.rs], show.rs)
    }
    if (show[3]) {
        sqrtabsr <- sqrt(abs(rs))
        ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
        yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
        yhn0 <- if (is.null(w)) 
            yh
        else yh[w != 0]
        plot(yhn0, sqrtabsr, xlab = l.fit, ylab = yl, main = main, 
            ylim = ylim, type = "n", ...)
        panel(yhn0, sqrtabsr, ...)
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(caption[3], 3, 0.25)
        if (id.n > 0) 
            text.id(yhn0[show.rs], sqrtabsr[show.rs], show.rs)
    }
    if (show[4]) {
        if (id.n > 0) {
            show.r <- order(-cook)[iid]
            ymx <- cook[show.r[1]] * 1.075
        }
        else ymx <- max(cook, na.rm = TRUE)
        plot(cook, type = "h", ylim = c(0, ymx), main = main, 
            xlab = "Obs. number", ylab = "Cook's distance", ...)
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(caption[4], 3, 0.25)
        if (id.n > 0) 
            text.id(show.r, cook[show.r], show.r, adj.x = FALSE)
    }
    if (show[5]) {
        ylim <- range(rs, na.rm = TRUE)
        if (id.n > 0) {
            ylim <- extendrange(r = ylim, f = 0.08)
            show.r <- order(-cook)[iid]
        }
        if (isConst.hat) {
            caption[5] <- "Constant Leverage:\n Residuals vs Factor Levels"
            aterms <- attributes(terms(x))
            dcl <- aterms$dataClasses[-aterms$response]
            facvars <- names(dcl)[dcl %in% c("factor", "ordered")]
            mf <- model.frame(x)[facvars]
            effM <- mf
            for (j in seq(length = ncol(mf))) effM[, j] <- sapply(split(yh, 
                mf[, j]), mean)[mf[, j]]
            dm <- data.matrix(mf)[do.call(order, effM), , drop = FALSE]
            nf <- length(nlev <- unlist(unname(lapply(x$xlevels, 
                length))))
            ff <- if (nf == 1) 
                1
            else rev(cumprod(c(1, nlev[nf:2])))
            xx <- facval <- (dm - 1) %*% ff
            plot(facval, rs, xlim = c(-1/2, sum((nlev - 1) * 
                ff) + 1/2), ylim = ylim, xaxt = "n", main = main, 
                xlab = "Factor Level Combinations", ylab = ylab23, 
                type = "n", ...)
            axis(1, at = ff[1] * (1:nlev[1] - 1/2) - 1/2, labels = x$xlevels[[1]][order(sapply(split(yh, 
                mf[, 1]), mean))])
            mtext(paste(facvars[1], ":"), side = 1, line = 0.25, 
                adj = -0.05)
            abline(v = ff[1] * (0:nlev[1]) - 1/2, col = "gray", 
                lty = "F4")
            panel(facval, rs, ...)
            abline(h = 0, lty = 3, col = "gray")
        }
        else {
            xx <- hatval
            xx[xx >= 1] <- NA
            plot(xx, rs, xlim = c(0, max(xx, na.rm = TRUE)), 
                ylim = ylim, main = main, xlab = "Leverage", 
                ylab = ylab23, type = "n", ...)
            panel(xx, rs, ...)
            abline(h = 0, v = 0, lty = 3, col = "gray")
            if (one.fig) 
                title(sub = sub.caption, ...)
            if (length(cook.levels)) {
                p <- length(coef(x))
                usr <- par("usr")
                hh <- seq(min(r.hat[1], r.hat[2]/100), usr[2], 
                  length = 101)
                for (crit in cook.levels) {
                  cl.h <- sqrt(crit * p * (1 - hh)/hh)
                  lines(hh, cl.h, lty = 2, col = 2)
                  lines(hh, -cl.h, lty = 2, col = 2)
                }
                legend("bottomleft", legend = "Cook's distance", 
                  lty = 2, col = 2, bty = "n")
                xmax <- min(0.99, usr[2])
                ymult <- sqrt(p * (1 - xmax)/xmax)
                aty <- c(-sqrt(rev(cook.levels)) * ymult, sqrt(cook.levels) * 
                  ymult)
                axis(4, at = aty, labels = paste(c(rev(cook.levels), 
                  cook.levels)), mgp = c(0.25, 0.25, 0), las = 2, 
                  tck = 0, cex.axis = cex.id, col.axis = 2)
            }
        }
        mtext(caption[5], 3, 0.25)
        if (id.n > 0) {
            y.id <- rs[show.r]
            y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
            text.id(xx[show.r], y.id, show.r)
        }
    }
    if (!one.fig && par("oma")[3] >= 1) 
        mtext(sub.caption, outer = TRUE, cex = 1.25)
    invisible()
}
