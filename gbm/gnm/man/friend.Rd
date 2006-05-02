\name{friend}
\alias{friend}
\docType{data}
\title{ Occupation of Respondents and Their Closest Friend }
\description{
  Cross-classification of the occupation of respondent and that of their
  closest friend. Data taken from wave 10 (year 2000) of the British
  Household Panel Survey.  
}
\usage{data(friend)}
\format{
 A table of counts, with classifying factors \code{r} (respondent's
 occupational category; levels \code{1:31}) and \code{c} (friend's
 occupational category; levels \code{1:31}).
}
\source{
  Chan, T.W. and Goldthorpe, J.H. (2004) Is there a status order 
  in contemporary British society: Evidence from the occupational 
  structure of friendship, \emph{European Sociological Review},
  \bold{20}, 383--401. 
}
\examples{
set.seed(1)
data(friend)

##  Fit an association model with homogeneous row-column effects
rc1 <- gnm(Freq ~ r + c + Diag(r,c) + Nonlin(MultHomog(r, c)),
           family = poisson, data = friend)
rc1

##  Extend to two-component interaction (using workaround)
original <- as.data.frame(friend)
r2 <- original$r
c2 <- original$c
rc2 <- update(rc1, . ~ . + Nonlin(MultHomog(r2, c2)))
rc2

}
\keyword{datasets}
