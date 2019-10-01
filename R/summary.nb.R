#  A generalized version of MASS::summary.negbin
#
#  Copyright (C) 1994-2014 W. N. Venables and B. D. Ripley
#  Copyright (C) 2019 Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
summary.nb <- function(object, dispersion = 1, correlation = FALSE, ...){
    if(is.null(dispersion)) dispersion <- 1
    summ <- c(summary.gnm(object, dispersion = dispersion,
                          correlation = correlation),
               object[c("theta", "SE.theta", "twologlik", "th.warn")])
    class(summ) <- c("summary.negbin", "summary.gnm")
    summ
}

vcov.nb <- function(object, dispersion = 1, with.eliminate = FALSE, ...){
    if(is.null(dispersion)) dispersion <- 1
    vcov.gnm(object, dispersion = dispersion, 
             with.eliminate = with.eliminate, ...)
}