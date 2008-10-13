\name{plot.gnm}
\alias{plot.gnm}
\title{ Plot Diagnostics for a gnm Object }
\description{
  Five plots are available: a plot of residuals against fitted values, a
  Scale-Location plot of sqrt{| residuals |} against fitted values, a
  Normal Q-Q plot, a plot of Cook's distances versus row labels, and a
  plot of residuals against leverages.  By default, all except the
  fourth are produced.
}
\usage{
\method{plot}{gnm}(x, which = c(1:3, 5), caption = c("Residuals vs Fitted",
         "Normal Q-Q", "Scale-Location", "Cook's distance",
         "Residuals vs Leverage"),
         panel = if (add.smooth) panel.smooth else points,
         sub.caption = NULL, main = "",
         ask = prod(par("mfcol")) < length(which) && dev.interactive(),
         ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
         qqline = TRUE, cook.levels = c(0.5, 1),
         add.smooth = getOption("add.smooth"), label.pos = c(4, 2),
         cex.caption = 1)    
}
\arguments{
  \item{x}{ a \code{"gnm"} object. }
  \item{which}{ a subset of the numbers 1:5 specifying which plots to
    produce (out of those listed in Description section). }
  \item{caption}{ captions to appear above the plots. }
  \item{panel}{ panel function.  The useful alternative to \code{points},
          \code{panel.smooth} can be chosen by \code{add.smooth = TRUE}. }
  \item{sub.caption}{ common title - above figures if there are
    multiple; used as \code{sub} (s.\code{title}) otherwise. If
    \code{NULL}, as by default, a possible shortened version of
    \code{deparse(x$call)} is used. }
  \item{main}{ title to each plot - in addition to the above \code{caption}. }
  \item{ask}{ logical; if \code{TRUE}, the user is asked before each
    plot, see \code{par(ask = .)}.}
  \item{\dots}{ other parameters to be passed through to plotting functions. }
  \item{id.n}{ number of points to be labelled in each plot starting
    with the most extreme. }
  \item{labels.id}{ vector of labels, from which the labels for extreme
    points will be chosen. \code{NULL} uses observation numbers. }
  \item{cex.id}{ magnification of point labels. }
  \item{qqline}{ logical indicating if a \code{qqline()} should be added
    to the normal Q-Q plot.}
  \item{cook.levels}{ levels of Cook's distance at which to draw contours. }
  \item{add.smooth}{ logical indicating if a smoother should be added to
    most plots; see also \code{panel} above.}
  \item{label.pos}{ positioning of labels, for the left half and right
    half of the graph respectively, for plots 1-3. }
  \item{cex.caption}{ controls the size of 'caption'. }
}
\details{
  \code{sub.caption} - by default the function call - is shown as a subtitle
  (under the x-axis title) on each plot when plots are on separate
  pages, or as a subtitle in the outer margin (if any) when there
  are multiple plots per page.

  The "Scale-Location" plot, also called "Spread-Location" or "S-L"
  plot, takes the square root of the absolute residuals in order to
  diminish skewness (sqrt{| E |} is much less skewed than | E | for
  Gaussian zero-mean E).

  The S-L, the Q-Q, and the Residual-Leverage plot, use
  \emph{standardized} residuals which have identical variance (under the
  hypothesis).  They are given as R[i] / (s*sqrt(1 - h.ii)) where
  h.ii are the diagonal entries of the hat matrix,
  \code{influence()$hat}, see also \code{\link{hat}}. 

  The Residual-Leverage plot shows contours of equal Cook's
  distance, for values of \code{cook.levels} (by default 0.5 and 1) and
  omits cases with leverage one. If the leverages are constant, as
  typically in a balanced \code{aov} situation, the plot uses factor
  level combinations instead of the leverages for the x-axis.
}
\author{ Heather Turner (adaptation of \code{plot.lm}, as of R 2.2.1) }
\seealso{ \code{\link{gnm}}, \code{\link{plot.lm}} }
\examples{
set.seed(1)
data(occupationalStatus)

##  Fit an association model with homogeneous row-column effects
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), family = poisson,
               data = occupationalStatus)

## Plot model diagnostics
plot(RChomog)

## Put 4 plots on 1 page; allow room for printing model formula in outer margin:
par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
title <- paste(deparse(RChomog$formula, width.cutoff = 50), collapse = "\n")
plot(RChomog, sub.caption = title)

## Fit smoother curves
plot(RChomog, sub.caption = title, panel = panel.smooth)
plot(RChomog, sub.caption = title, panel = function(x,y) panel.smooth(x, y, span = 1))

}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ hplot }
