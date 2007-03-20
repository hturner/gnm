##################################################################
FLog <- function(x) {
## This always returns a finite number.
## Used to stop iteration routine from failing.
log( pmax(1e-100, x) )
}
##################################################################
Sclgstc <- function(x) {
## Used to ensure transformed parameter is non-negative
## and finite, so iteration routine doesn’t fail.
1e5/(1+exp(-x))
}
##################################################################

##################################################################
andyBeli <- function(x, Fpeakage=~1, Fpeak.ht=~1,
                 Fpksharp=~1, Fleft.ep=~1) {
    ## See coding for Bell. Here we have
    ## h(x) = i(peakage) * log( i(peakage)/i(x) ) + x - peakage
    ## Note that this is a reparameterisation of
    ## c + AL*log(x-start.point) - AR*x
    ## This is equivalent to Bell with end.point = +Inf
    ## get design matrices for dependencies
    gnmData <- getModelFrame()
    Mpeakage <- model.matrix(Fpeakage, data=gnmData)
    Mpeak.ht <- model.matrix(Fpeak.ht, data=gnmData)
    Mpksharp <- model.matrix(Fpksharp, data=gnmData)
    Mleft.ep <- model.matrix(Fleft.ep, data=gnmData)
    ## create index and labels for parameters
    Npeakage <- ncol(Mpeakage)
    Npeak.ht <- ncol(Mpeak.ht)
    Npksharp <- ncol(Mpksharp)
    Nleft.ep <- ncol(Mleft.ep)
    id <- rep(1:4, c(Npeakage, Npeak.ht, Npksharp, Nleft.ep))
    if (Npeakage>1)
        Lpeakage <- paste("peakage", colnames(Mpeakage), sep=".")
    else
        Lpeakage <- "peakage"
    if (Npeak.ht>1)
        Lpeak.ht <- paste("peak.ht", colnames(Mpeak.ht), sep=".")
    else
        Lpeak.ht <- "peak.ht"
    if (Npksharp>1)
        Lpksharp <- paste("pksharp", colnames(Mpksharp), sep=".")
    else
        Lpksharp <- "pksharp"
    if (Nleft.ep>1)
        Lleft.ep <- paste("left.ep", colnames(Mleft.ep), sep=".")
    else
        Lleft.ep <- "left.ep"
    labels <- c(Lpeakage, Lpeak.ht, Lpksharp, Lleft.ep)
    env <- environment()
    ## function to calculate linear predictor
    ## (and some intermediate results for partial derivatives)
    predictor <- function(coef) {
        assign("peakage", drop(Mpeakage %*% coef[id==1]), env)
        assign("peak.ht", drop(Mpeak.ht %*% coef[id==2]), env)
        assign("pksharp", drop(Mpksharp %*% coef[id==3]), env)
        assign("left.ep", drop(Mleft.ep %*% coef[id==4]), env)
        assign("S",peakage-min(x)+1e-5+Sclgstc(left.ep), env)
        assign("X", S + x - peakage, env)
        assign("H", S - 5, env)
        assign("SX", S/X, env)
        assign("SH", S/H, env)
        assign("lSX", FLog(SX), env)
        assign("lSH", FLog(SH), env)
        assign("v", S * lSH - 5, env)
        assign("V", (S * lSX - S + X)/v, env)
        assign("p", peak.ht-exp(pksharp)*V, env)
        p
    }
    ## function to calculate partial derivatives
    localDesignFunction <- function(coef, predictor, ...) {
        cbind(Mpeakage * exp(pksharp) / v *
              (V*(1-SH+lSH)-lSX),
              Mpeak.ht,
              Mpksharp * (p-peak.ht),
              Mleft.ep * Sclgstc(left.ep) / (1+exp(left.ep)) *
              exp(pksharp) / v *
              (V*(1-SH+lSH)-1+SX-lSX))
    }
    list(start = rep(c(25,-2.25,-1.5,-12), c(Npeakage,
         Npeak.ht, Npksharp, Nleft.ep) ),
         labels = labels, predictor = predictor,
         localDesignFunction = localDesignFunction)
}
##################################################################
andyBeliVariables <- function(x, Fpeakage=~1, Fpeak.ht=~1,
                          Fpksharp=~1, Fleft.ep=~1) {
    .call <- match.call()
    sapply(c(.call$x, Fpeakage[[2]], Fpeak.ht[[2]], Fpksharp[[2]],
           Fleft.ep[[2]]), deparse)
}
##################################################################
