\name{occupationalStatus}
\alias{occupationalStatus}
\docType{data}
\title{ Data on Occupational Status of Fathers and their Sons }
\description{
  Cross-classification of a sample of British males according to each
  subject's occupational status and his father's occupational status.
}

\usage{data(occupationalStatus)}
\format{
  A table of counts, with classifying factors \code{origin} (father's
  occupational status; levels \code{1:8}) and \code{destination} (son's
  occupational status; levels \code{1:8}).
}
\source{
  Goodman, L. A. (1979) Simple Models for the Analysis of Association in
  Cross-Classifications having Ordered
  Categories. \emph{J. Am. Stat. Assoc.}, \bold{74(367)}, 537--552.
}
\examples{
set.seed(1)
data(occupationalStatus)

##  Fit a uniform association model separating diagonal effects
Rscore <- scale(as.numeric(row(occupationalStatus)), scale = FALSE)
Cscore <- scale(as.numeric(col(occupationalStatus)), scale = FALSE)
Uniform <- glm(Freq ~ origin + destination + Diag(origin, destination) + 
               Rscore:Cscore, family = poisson, data = occupationalStatus)

##  Fit an association model with homogeneous row-column effects
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               Nonlin(MultHomog(origin, destination)), family = poisson,
               data = occupationalStatus)

##  Fit an association model with separate row and column effects
RC <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
          Mult(origin, destination), family = poisson,
          data = occupationalStatus)

##  Partition row and column effects
row.and.column.effects <- c(df = Uniform$df.res - RC$df.res,
                            deviance = Uniform$dev - RC$dev)
homogeneous.row.column.effect <- c(df = Uniform$df.res - RChomog$df.res,
                                    deviance = Uniform$dev - RChomog$dev)
row.column.difference.effect <- row.and.column.effects -
    homogeneous.row.column.effect
rbind(homogeneous.row.column.effect,
      row.column.difference.effect,
      row.and.column.effects)
}
\keyword{datasets}
