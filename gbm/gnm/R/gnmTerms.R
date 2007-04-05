gnmTerms <- function(formula, eliminate = NULL, data = NULL)
{
    if (!is.null(substitute(e, list(e = eliminate)))) {
        tmp <- .Internal(update.formula(formula,
                                        substitute(~ -1 + e + .,
                                                   list(e = eliminate))))
        environment(tmp) <- environment(formula)
        data <- data[!names(data) %in% deparse(eliminate)]
        formula <- formula(terms.formula(tmp, simplify = TRUE,
                                         keep.order = TRUE, data = data))
    }
    fullTerms <- terms(formula, specials = c("Const", "instances"),
                       simplify = TRUE, keep.order = TRUE, data = data)

    if (is.empty.model(fullTerms))
        return(fullTerms)

    inst <- attr(fullTerms, "specials")$instances
    if (length(inst)) {
        termLabels <- c("0"[!attr(fullTerms, "intercept")],
                        attr(fullTerms, "term.labels"))
        instLabels <- as.list(attr(fullTerms, "variables"))[inst + 1]
        termLabels[termLabels %in% instLabels] <- sapply(instLabels, eval)
        variables <- as.character(attr(fullTerms, "variables"))[-1]
        offsetLabels <- variables[attr(fullTerms, "offset")]
        response <- variables[attr(fullTerms, "response")][1][[1]]
        fullTerms <- terms(reformulate(c(termLabels, offsetLabels), response),
                           keep.order = TRUE, data = data, specials = "Const")
    }

    termLabels <- c("1"[attr(fullTerms, "intercept")],
                    attr(fullTerms, "term.labels"))
    variables <- predvars <- as.list(attr(fullTerms, "variables"))[-1]
    specials <- which(sapply(variables, function(x) {
        length(x) > 1 && inherits(match.fun(x[[1]]), "nonlin")
    }))
    if (!length(specials)) {
        if (is.null(eliminate))
            return(fullTerms)
        n <- length(termLabels)
        attributes(fullTerms) <-
            c(attributes(fullTerms),
              list(unitLabels = termLabels,
                   common = logical(n),
                   block = numeric(n),
                   match = !logical(n),
                   assign = rep(1, n),
                   classID = rep.int("Linear", n),
                   NonlinID = character(n),
                   prefixLabels = character(n),
                   varLabels = termLabels,
                   predictor = lapply(termLabels, as.name),
                   class = c("gnmTerms", "terms", "formula")))
        return(fullTerms)
    }

    specialTerms <- rownames(attr(fullTerms, "factors"))[specials]
    specialTerms <- strsplit(specialTerms, ", inst = |,? ?\\)$", perl = TRUE)
    term <- sapply(specialTerms, "[", 1)
    inst <- as.numeric(sapply(specialTerms, "[", 2))

    patch <- term %in% term[inst > 1] & is.na(inst)
    termLabels[termLabels %in% specials[patch]] <-
        paste(term[patch], ", inst = 1)")
    inst[patch] <- 1

    nonsense <- tapply(inst, term, FUN = function(x)
                   {!is.na(x) && !identical(as.integer(x), seq(x))})
    if (any(nonsense))
        stop("Specified instances of ",
             paste(names(nonsense)[nonsense], ")"),
             " are not in sequence")

    const <- attr(fullTerms, "specials")$Const
    if (length(const)) {
        termLabels <- termLabels[!termLabels %in% variables[const]]
        predvars[const] <- lapply(variables[const], eval)
    }
    offsetVars <- variables[c(attr(fullTerms, "offset"), const)]
    nonlinear <- termLabels %in% variables[specials]
    variables <- variables[-specials]
    predvars <- predvars[-specials]

    unitLabels <- varLabels <- as.list(termLabels)
    predictor <- lapply(termLabels, as.name)
    names(predictor) <- unitLabels
    n <- length(unitLabels)
    blockList <- as.list(numeric(n))
    match <- as.list(!logical(n))
    common <- as.list(logical(n))
    class <- as.list(rep.int("Linear", n))
    NonlinID <- prefixLabels <- as.list(character(n))
    start <- vector("list", n)
    adj <- 1

    for (j in which(nonlinear)) {
        if (identical(substr(unitLabels[[j]], 0, 7), "Nonlin(")){
            tmp <- eval(parse(text = unitLabels[[j]]))
            tmp$varLabels <- tmp$unitLabels <- unitLabels[[j]]
            tmp$predictor <- paste("`", unitLabels[[j]], "`", sep = "")
        }
        else
            tmp <- do.call("nonlinTerms",
                           eval(parse(text = unitLabels[[j]]), data))
        unitLabels[[j]] <- tmp$unitLabels
        if (!identical(tmp$prefix, "#")) {
            bits <- hashSplit(tmp$prefix)
            if (length(bits) > 1) {
                n <- length(tmp$hashLabels)
                matched <- tmp$matchID > 0 & !duplicated(tmp$matchID)
                dot <- (tmp$hashLabels[matched])[order(tmp$matchID[matched])]
                prefix <- matrix(dot, max(tmp$matchID), n)
                prefix[cbind(tmp$matchID, seq(n))] <- "."
                prefix <- rbind(character(n), prefix)
                sep <- rep(".", n)
                sep[!tmp$matchID] <- ""
                prefixLabels[[j]] <- paste(apply(prefix, 2, paste, bits,
                                                 sep = "", collapse = ""),
                                           sep, tmp$suffix, sep = "")
                for (i in unique(tmp$common[duplicated(tmp$common)])) {
                    dotCommon <- dot
                    commonID <- tmp$common == i
                    dotCommon[tmp$matchID[commonID]] <- "."
                    prefixLabels[[j]][commonID] <-
                        paste(paste(c("", dotCommon), bits, sep = "",
                                    collapse = ""),
                              tmp$suffix[commonID], sep[commonID],
                              paste(tmp$unitLabels[commonID], collapse = "|"),
                              sep = "")
                }
            }
            else
                prefixLabels[[j]] <- paste(tmp$prefix, tmp$suffix, sep = "")
        }
        else
            prefixLabels[[j]] <- tmp$varLabels

        varLabels[[j]] <- gsub("#", j, tmp$varLabels)
        predictor[[j]] <- parse(text = gsub("#", j, tmp$predictor))[[1]]
        blockList[[j]] <- tmp$block + adj
        match[[j]] <- as.logical(tmp$matchID)
        common[[j]] <- tmp$common %in% tmp$common[duplicated(tmp$common)]
        class[[j]] <- tmp$classID
        NonlinID[[j]] <- tmp$NonlinID
        start[j] <- list(tmp$start)
        adj <- max(c(0, blockList[[j]])) + 1
        variables <- c(variables, tmp$variables)
        predvars <- c(predvars, tmp$predvars)
    }

    if (length(predvars) > 1)
        nObs <- call("length", predvars[[1]])
    else if (!is.null(data))
        nObs <- call("length", as.name(names(data)[1]))
    else
        nObs <- 1

    attributes(fullTerms) <-
        c(attributes(fullTerms),
          list(offset = which(unique(variables) %in% offsetVars),
               variables = as.call(c(quote(list), unique(variables))),
               predvars = {do.call("substitute",
                                   list(as.call(c(quote(list),
                                                  unique(predvars))),
                                        list(nObs = nObs)))},
               unitLabels = unlist(unitLabels),
               common = unlist(common),
               block = unlist(blockList),
               match = unlist(match),
               assign = rep(seq(class), sapply(class, length)),
               classID = unlist(class),
               NonlinID = unlist(NonlinID),
               prefixLabels = unlist(prefixLabels),
               varLabels = unlist(varLabels),
               start = start,
               predictor = predictor,
               class = c("gnmTerms", "terms", "formula")))
    fullTerms
}
