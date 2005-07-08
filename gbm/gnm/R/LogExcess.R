LogExcess <- function(..., side = "left", constraint = NULL, method = "default") { # would be able to put constraint on x if specified
    gnmData <- getModelFrame()
    
    x <- gnmData[, as.character(substitute(...))]
    constraint <- with(gnmData, constraint)

    if (!side %in% c("left", "right")) {
        warning("'side' not \"right\", so assuming \"left\"", call. = FALSE)
        side <- "left"
    }
    switch(side,
           "left" =  if (constraint >= min(x)) {
               warning("'constraint' invalid, replacing with ", min(x) - 1,
                       "\n", call. = FALSE)
               constraint <- min(x) - 1
           },
           "right" = if (constraint <= max(x)) {
               warning("'constraint' invalid, replacing with ", max(x) + 1,
                       "\n", call. = FALSE)
               constraint <- max(x) + 1
           })
    
    x <- switch(side, "left" = x - constraint, "right" = constraint - x)

    if (method == "default")
        start <- NULL
    else if (method == "offsetNone")
        start <- function(n, family, offset) {
            c(suppressWarnings(glm.fit(log(x + 1), model.response(gnmData),
                                       weights = model.weights(gnmData),
                                       offset = model.offset(gnmData),
                                       family = family))$coef,
              0)
        }
    else if (method == "offsetOthers")
        start <- function(n, family, offset){
            if (!is.null(model.offset(gnmData)))
                offset <- offset + model.offset(gnmData)
            c(suppressWarnings(glm.fit(log(x + 1), model.response(gnmData),
                                       weights = model.weights(gnmData),
                                       offset = offset, family = family))$coef,
              0)
        }
    
    predictor <- function(coef) {
        coef[1]*log(exp(coef[2]) + x)
    }

    localDesignFunction <- function(coef, ...) {
        cbind(log(x + exp(coef[2])), coef[1] * exp(coef[2])/(x + exp(coef[2])))
    }

    list(start = start, labels = c("slope", "endpoint"), predictor = predictor,
         localDesignFunction = localDesignFunction)
}
