bell <- function(..., lowerMax = min(x) - 1, upperMin = max(x) + 1,
                 model = "full", b1 = NULL, b2 = NULL) {
    gnmData <- getModelFrame()
    x <- gnmData[, as.character(substitute(...))]

    ## adjust x
    left <- x - with(gnmData, lowerMax)
    right <- with(gnmData, upperMin) - x

    predictor <- function(coef) {
        switch(model,
               "full" = matrix(c(coef[1]*log(exp(coef[2]) + left),
               coef[3]*log(exp(coef[4]) + right)), ncol = 2),
               "fixedEndpoints" = matrix(c(coef[1]*log(left),
               coef[2]*log(right)), ncol = 2),
               "fixedSlopes" = matrix(c(b1 * log(exp(coef[1]) + left),
               b2 * log(exp(coef[2]) + right)), ncol = 2),
               "firstTerm" = matrix(coef[1]*log(exp(coef[2]) + left), ncol = 1))
    }

    localDesignFunction <- function(coef, predictor) {
        #slopeCol <- sweep(predictor, 2, coef[c(1, 3)], "/")
        #adjustCol <- sweep(1/exp(slopeCol), 2,
                           #coef[c(1, 3)] * exp(coef[c(2, 4)]), "*")
        #cbind(slopeCol[,1], adjustCol[,1], slopeCol[,2], adjustCol[,2])
        dLS <- switch(model,
                      "full" = log(left + exp(coef[2])),
                      "fixedEndpoints" = log(left),
                      "firstTerm" = log(left + exp(coef[2])),
                      NULL)
        dLA <- switch(model,
                      "full" = coef[1] * exp(coef[2])/(left + exp(coef[2])),
                      "fixedSlopes" = b1 * exp(coef[1])/(left + exp(coef[1])),
                      "firstTerm" = coef[1] * exp(coef[2])/(left + exp(coef[2])),
                      NULL)
        dRS <- switch(model,
                      "full" = log(right + exp(coef[4])),
                      "fixedEndpoints" = log(right),
                      NULL) 
        dRA <- switch(model,
                      "full" = coef[3] * exp(coef[4])/(right + exp(coef[4])),
                      "fixedSlopes" = b2 * exp(coef[2])/(right + exp(coef[2])),
                      NULL)
        cbind(dLS, dLA, dRS, dRA)
    }

    labels <- switch(model,
                     "full" = c("leftSlope", "leftAdjust", "rightSlope", "rightAdjust"),
                     "fixedEndpoints" = c("leftSlope", "rightSlope"),
                     "fixedSlopes" = c("leftAdjust", "rightAdjust"),
                     "firstTerm" = c("leftSlope", "leftAdjust"))
    list(labels = labels,
         predictor = predictor, localDesignFunction = localDesignFunction)
}
