\name{occupational.status}
\alias{occupational.status}
\docType{data}
\title{ Occupational Status of Fathers and their Sons }
\description{
  Cross-classification of a sample of British males according to each
  subject's occupational status and his father's occupational status.
}

\usage{data(occupational.status)}
\format{
  A table of counts, with classifying factors \code{origin} (father's
  occupational status; levels \code{1:8}) and \code{destination} (son's
  occupational status; levels \code{1:8}).
}
\source{
Goodman, L. A. (1979) Simple Models for the Analysis of Association in
Cross-Classifications having Ordered
Categories. \emph{J. Am. Stat. Assoc.}, \bold{74(367)}, 537-552.
}
\examples{
data(occupational.status)

##  Fit a uniform association model separating diagonal effects
Rscore <- scale(as.numeric(row(occupational.status)), scale = FALSE)
Cscore <- scale(as.numeric(col(occupational.status)), scale = FALSE)
Uniform <- glm(Freq ~ origin + destination + Diag(origin, destination) + 
  Rscore*Cscore, family = poisson, data = occupational.status)

##  Fit an association model with homogeneous row-column effects
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
  Nonlin(MultHomog(origin, destination)), family = poisson,
  data = occupational.status)

##  Fit an association model with separate row and column effects
RC <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
  Mult(origin, destination), family = poisson, data = occupational.status,
  start = c(coef(RChomog),
      coef(RChomog)[(length(coef(RChomog)) - 7):length(coef(RChomog))]))

##  Partition row and column effects
row.and.column.effects <- c(df = Uniform$df.res - RC$df.res,
                            deviance = Uniform$dev - RC$dev)
homogeneous.row.column.effects <- c(df = Uniform$df.res - RChomog$df.res,
                                    deviance = Uniform$dev - RChomog$dev)
row.column.effects <- row.and.column.effects -
                      homogeneous.row.column.effects
rbind(homogeneous.row.column.effects,
      row.column.effects,
      row.and.column.effects)
}
\keyword{datasets}
