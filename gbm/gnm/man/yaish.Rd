\name{yaish}
\alias{yaish}
\docType{data}
\title{ Class mobility by education in Israel}
\description{
  A 3-way contingency table of father/son pairs, classified by father's
  social class (\code{orig}), son's social class (\code{dest}) and son's 
  education level (\code{educ}).
}
\usage{data(yaish)}
\format{
A table of counts, with classifying factors \code{educ} (levels \code{1:5}),
\code{orig} (levels \code{1:7}) and \code{dest} (levels \code{1:7}).
}
\references{
Yaish, M (1998). Opportunities, Little Change. Class Mobility in Israeli Society: 
1974-1991.  D.Phil. Thesis, Nuffield College, University of Oxford.
}
\examples{
data(yaish)
yaish.frame <- as.data.frame.table(yaish)
## Fit the "UNIDIFF" mobility model across education levels
yaish.unidiff <- gnm(Freq ~ educ:orig + educ:dest + Mult(Exp(educ), orig:dest), family = poisson, data = yaish)
## Look at a table of model residuals
residuals(yaish.unidiff)
}
\keyword{datasets}
