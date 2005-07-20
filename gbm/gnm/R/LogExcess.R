LogExcess <- function(x, side = "left", constraint = ifelse(side == "right",
                                        max(x) + 1, min(x) - 1)) {
    
    if (!side %in% c("left", "right")) {
        warning("'side' not \"right\", so assuming \"left\"")
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
    
    predictor <- function(coef) {
        coef[1]*log(exp(coef[2]) + x)
    }

    localDesignFunction <- function(coef, ...) {
        cbind(log(x + exp(coef[2])), coef[1] * exp(coef[2])/(x + exp(coef[2])))
    }

    list(start = c(NA, 0), labels = c("slope", "excess"),
         predictor = predictor, localDesignFunction = localDesignFunction)
}
