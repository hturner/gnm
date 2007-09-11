DrefWeights <- function(model) {
    ind <- pickCoef(model, "delta[1-9]")
    if (any(!checkEstimable(model, diag(seq(along = coef(model)))[,ind]),
        na.rm = TRUE)){
        message("Refitting with parameters of first Dref weight constrained ",
                "to zero")
        constrain <- pickCoef(model, "delta1")
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
        X <- matrix(1)
    }
    nw <- length(preds)
    nmod <- nrow(X)
    delta <- matrix(parameters(model)[ind], nmod)
    ind <- c(t(matrix(ind, nmod, nw)))
    vcovDelta <- vcov(model)[ind, ind, drop = FALSE]
    wc <- 1/rowSums(exp(X %*% delta))
    wu <- exp(X %*% delta)*wc
    XX <- matrix(apply(X, 2, rep, nw), nrow(X))
    out <- list()
    for (i in 1:nw) {
        d <- -wu[,i] * wu
        d <- c(wu[,i] * (col(wu) == i) + d) * XX
        se <- sqrt(rowSums(d %*% vcovDelta * d))
        out[[i]] <- drop(cbind(dat, weight = wu[,i], se = se))
    }
    names(out) <- as.character(preds)
    out
}
