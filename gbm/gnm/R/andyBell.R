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
andyBell <- function(x, Fpeakage=~1, Fpeak.ht=~1, Fpksharp=~1,
                 Fleft.ep=~1, Frightep=~1) {
    ## This is a Custom Plug-in Function written to enable gnm to
    ## fit nonlinear terms. See 3) Nonlinear terms in Generalized
    ## nonlinear models in R: an overview of the gnm package.
    ## This coding, adapted from coding by H Turner, deals with
    ## terms of the form peak.ht - exp(pksharp) * V where
    ## V = h(x) / h(peakage-5),
    ## h(x) = i(peakage) * log( i(peakage)/i(x) )
    ## +j(peakage) * log( j(peakage)/j(x) )
    ## i(x) = x - start.point
    ## j(x) = end.point - x
    ## start.point = min(x) - 1e-5 - Sclgstc(left.ep)
    ## end.point = max(x)+ 1e-5 + Sclgstc(rightep)
    ## Note that this is a reparameterisation of
    ## c + AL*log(x-start.point) + AR*log(end.point-x)
    ## get design matrices for dependencies
    gnmData <- getModelFrame()
    Mpeakage <- model.matrix(Fpeakage, data=gnmData)
    Mpeak.ht <- model.matrix(Fpeak.ht, data=gnmData)
    Mpksharp <- model.matrix(Fpksharp, data=gnmData)
    Mleft.ep <- model.matrix(Fleft.ep, data=gnmData)
    Mrightep <- model.matrix(Frightep, data=gnmData)
    ## create index and labels for parameters
    Npeakage <- ncol(Mpeakage)
    Npeak.ht <- ncol(Mpeak.ht)
    Npksharp <- ncol(Mpksharp)
    Nleft.ep <- ncol(Mleft.ep)
    Nrightep <- ncol(Mrightep)
    id <- rep(1:5, c(Npeakage, Npeak.ht, Npksharp,
                     Nleft.ep, Nrightep))
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
    if (Nrightep>1)
        Lrightep <- paste("rightep", colnames(Mrightep), sep=".")
    else
        Lrightep <- "rightep"
    labels <- c(Lpeakage, Lpeak.ht, Lpksharp, Lleft.ep, Lrightep)
    env <- environment()
    ## function to calculate linear predictor
    ## (and some intermediate results for partial derivatives)
    predictor <- function(coef) {
        assign("peakage", drop(Mpeakage %*% coef[id==1]), env)
        assign("peak.ht", drop(Mpeak.ht %*% coef[id==2]), env)
        assign("pksharp", drop(Mpksharp %*% coef[id==3]), env)
        assign("left.ep", drop(Mleft.ep %*% coef[id==4]), env)
        assign("rightep", drop(Mrightep %*% coef[id==5]), env)
        assign("S",
               list(L=peakage-min(x)+1e-5+Sclgstc(left.ep),
                    R=max(x)-peakage+1e-5+Sclgstc(rightep)), env)
        assign("X", list(L=S$L+x-peakage, R=S$R+peakage-x), env)
        assign("H", list(L=S$L-5, R=S$R+5), env)
        assign("SX", list(L=S$L/X$L, R=S$R/X$R), env)
        assign("SH", list(L=S$L/H$L, R=S$R/H$R), env)
        assign("lSX", list(L=FLog(SX$L), R=FLog(SX$R)), env)
        assign("lSH", list(L=FLog(SH$L), R=FLog(SH$R)), env)
        assign("v", S$L*lSH$L+S$R*lSH$R, env)
        assign("V", (S$L*lSX$L+S$R*lSX$R)/v, env)
        assign("p", peak.ht-exp(pksharp)*V, env)
        p
    }
    ## function to calculate partial derivatives
    localDesignFunction <- function(coef, predictor, ...) {
        cbind(Mpeakage * exp(pksharp) / v *
              (V*(SH$R-SH$L+lSH$L-lSH$R)+lSX$R-lSX$L),
              Mpeak.ht,
              Mpksharp * (p-peak.ht),
              Mleft.ep * Sclgstc(left.ep) / (1+exp(left.ep)) *
              exp(pksharp) / v *
              (V*(lSH$L+1-SH$L)-1+SX$L-lSX$L),
              Mrightep * Sclgstc(rightep) / (1+exp(rightep)) *
              exp(pksharp) / v *
              (V*(lSH$R+1-SH$R)-1+SX$R-lSX$R) )
    }
    list(start = rep(c(25,-2.25,-1.5,-12,-10), c(Npeakage,
         Npeak.ht, Npksharp, Nleft.ep, Nrightep) ),
         labels = labels, predictor = predictor,
         localDesignFunction = localDesignFunction)
}
##################################################################
andyBellVariables <- function(x, Fpeakage=~1, Fpeak.ht=~1, Fpksharp=~1,
                          Fleft.ep=~1, Frightep=~1) {
    .call <- match.call()
    sapply(c(.call$x, Fpeakage[[2]], Fpeak.ht[[2]], Fpksharp[[2]],
           Fleft.ep[[2]], Frightep[[2]]), deparse)
}
##################################################################
