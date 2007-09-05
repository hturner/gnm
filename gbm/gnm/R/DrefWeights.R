DrefWeights <- function(model) {
    constrain <- pickCoef(model, "delta1")
    ind <- setdiff(pickCoef(model, "delta[2-9]"), constrain)
    if (any(!checkEstimable(model, diag(seq(along = coef(model)))[,ind]))) {
        message("Refitting with first Dref weight constrained to zero")
        model <- update(model, constrain = constrain, start = coef(model),
                        trace = FALSE, verbose = FALSE)
    }
    t <- terms(formula(model), specials = "Dref")
    DrefCall <- attr(t, "variables")[[attr(t, "specials")$Dref + 1]]
    preds <- match.call(Dref, DrefCall, expand.dots = FALSE)[["..."]]
    formula <- as.formula(DrefCall$delta)
    if (length(formula)) {
        dat <- model.frame(formula, data = model.frame(model))
        X <- unique(model.matrix(formula, data = dat))
        dat <- dat[rownames(X), , drop = FALSE]
        rownames(dat) <- rownames(X) <- NULL
    }
    else {
        dat <- numeric(0)
        X <- matrix(1, 1, length(ind))
    }
    nw <- length(preds)
    nmod <- nrow(X)
    delta <- matrix(coef(model)[ind], nmod)
    ind <- c(t(matrix(ind, nmod, nw - 1)))
    vcovDelta <- vcov(model)[ind, ind, drop = FALSE]
    wc <- unname(1 + rowSums(exp(X %*% delta)))^(-1)
    wu <- exp(X %*% delta)*wc
    XX <- matrix(apply(X, 2, rep, nw - 1), nrow(X))
    out <- list()
    d <- -wc * (1 - wc) * XX
    se <- sqrt(rowSums(d %*% vcovDelta * d))
    out[[1]] <- cbind(dat, weight = wc, se = se)
    for (i in 1:(nw - 1)) {
        d <- -wu[,i] * wu
        d <- (d[,i] + wu[,i]) * XX
        se <- sqrt(rowSums(d %*% vcovDelta * d))
        out[[i + 1]] <- cbind(dat, weight = wu[,i], se = se)
    }
    names(out) <- as.character(preds)
    out
}
