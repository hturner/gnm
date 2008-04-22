"wedderburn" <-
    function (link = "logit") 
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link") 
            linktemp <- eval(link)
    }
    if (any(linktemp == c("logit", "probit", "cloglog"))) 
        stats <- make.link(linktemp)
    else stop(paste(linktemp,
                    "link not available for wedderburn quasi-family;", 
                    "available links are",
                    "\"logit\", \"probit\" and \"cloglog\""))
    variance <- function(mu)  mu^2 * (1-mu)^2
    validmu <- function(mu) {
        all(mu > 0) && all(mu < 1)}
    dev.resids <- function(y, mu, wt){
        eps <-  0.0005
        2 * wt * (y/mu + (1 - y)/(1 - mu) - 2 +
                  (2 * y - 1) * log((y + eps)*(1 - mu)/((1- y + eps) * mu)))
    }
    aic <- function(y, n, mu, wt, dev) NA
    initialize <- expression({
        if (any(y < 0 | y > 1)) stop(paste(
                   "Values for the wedderburn family must be in [0,1]"))
        n <- rep.int(1, nobs)
        mustart <- (y + 0.1)/1.2
    })
    structure(list(family = "wedderburn",
                   link = linktemp, 
                   linkfun = stats$linkfun,
                   linkinv = stats$linkinv,
                   variance = variance, 
                   dev.resids = dev.resids,
                   aic = aic,
                   mu.eta = stats$mu.eta, 
                   initialize = initialize,
                   validmu = validmu,
                   valideta = stats$valideta), 
              class = "family")
}
